#include <functional>
#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include "pgenlib_read.h"
#include "openmx.h"
#include "LoadDataAPI.h"

using namespace plink2;

struct BgenXfer {
	dataPtr dp;
	std::vector<double> prob;
	int dr;
	double total;
	int nrows;
	std::function<bool(int)> skipFn;
	BgenXfer(dataPtr &_dp, std::function<bool(int)> _skipFn) : dp(_dp), total(0), nrows(0), skipFn(_skipFn) {};
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {}
	void set_min_max_ploidy(genfile::bgen::uint32_t min_ploidy, genfile::bgen::uint32_t max_ploidy,
				genfile::bgen::uint32_t min_entries, genfile::bgen::uint32_t max_entries)
	{
		dr = 0;
		if (min_ploidy != 2 || max_ploidy != 2 || min_entries != 3 || max_entries != 3) {
			mxThrow("set_min_max_ploidy %u %u %u %u, not implemented",
				min_ploidy, max_ploidy, min_entries, max_entries);
		}
	}
	bool set_sample( std::size_t i ) { return !skipFn(i); }
	void set_number_of_entries(std::size_t ploidy,
				   std::size_t number_of_entries,
				   genfile::OrderType order_type,
				   genfile::ValueType value_type)
	{
		if( value_type != genfile::eProbability ) {
			mxThrow("value_type != genfile::eProbability");
		}
		prob.resize(number_of_entries);
	}
	void set_value( genfile::bgen::uint32_t entry_i, double value ) {
		prob[entry_i] = value;
		if (entry_i == 2) {
			double dosage = prob[1] * 1 + prob[0] * 2;
			if (dosage == 0.0) {
				dosage = NA_REAL;
			} else {
				total += dosage;
				nrows += 1;
			}
			dp.realData[dr++] = dosage;
		}
	}

	void set_value( genfile::bgen::uint32_t entry_i, genfile::MissingValue) {
		// Seems that NA is encoded as an exact zero?
		throw std::runtime_error("not implemented");
	}
};

struct LoadDataBGENProvider2 : public LoadDataProvider2<LoadDataBGENProvider2> {
	int cpIndex;
	genfile::bgen::View::UniquePtr bgenView;

	virtual const char *getName() { return "bgen"; };
	virtual void init(SEXP rObj1) {
		RObject rObj(rObj1);
		bool byrow = as<bool>(rObj.slot("byrow"));
		if (!byrow) mxThrow("byrow=FALSE is not implemented for bgen format");
		requireFile(rObj1);
	}
	virtual void loadRowImpl(int index);
	virtual void addCheckpointColumns(std::vector< std::string > &cp)
	{
		cpIndex = cp.size();
		cp.push_back("SNP");
		cp.push_back("RSID");
		cp.push_back("CHR");
		cp.push_back("BP");
		cp.push_back("A1");
		cp.push_back("A2");
		cp.push_back("MAF");
	}
	virtual int getNumVariants();
};

int LoadDataBGENProvider2::getNumVariants()
{
	if (bgenView.get() == 0) return 0;
	return bgenView->number_of_variants();
}

void LoadDataBGENProvider2::loadRowImpl(int index)
{
	// discard m_postheader_data? TODO
	if (columns.size() != 1) mxThrow("%s: bgen only has 1 column, not %d",
					 name, int(columns.size()));
	if (colTypes[0] != COLUMNDATA_NUMERIC) mxThrow("%s: bgen contains a numeric dosage", name);

	if (curRecord != index) bgenView.reset();

	if (bgenView.get() == 0) {
		using namespace genfile::bgen ;
		//using namespace Rcpp ;
		std::string bgen(filePath);
		std::string bgenIndex = bgen + ".bgi";
		bgenView = View::create( filePath ) ;
		auto query = IndexQuery::create( bgenIndex ) ;
		query->from_row(index);
		query->initialise();
		bgenView->set_query( query ) ;
		curRecord = index;
		if (srcRows != int(bgenView->number_of_samples())) {
			mxThrow("%s: %s has %d rows but %s has %d samples",
				name, dataName, srcRows, filePath.c_str(),
				int(bgenView->number_of_samples()));
		}
		loadCounter += 1;
	}

	std::string SNPID, rsid, chromosome ;
	genfile::bgen::uint32_t position ;
	std::vector< std::string > alleles ;
	if (!bgenView->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles )) {
		mxThrow("%s: %s has no more varients", name, filePath.c_str());
	}
	BgenXfer xfer(stripeData[0], [&](int rx)->bool{ return skipRow(rx); });
	bgenView->read_genotype_data_block(xfer);
	curRecord += 1;

	auto &rc = *rawCols;
	for (int cx=0; cx < int(columns.size()); ++cx) {
		rc[ columns[cx] ].setBorrow(stripeData[cx]);
	}
	if (checkpoint) {
		auto &cv = *checkpointValues;
		cv[cpIndex] = SNPID;
		cv[cpIndex+1] = rsid;
		cv[cpIndex+2] = chromosome;
		cv[cpIndex+3] = string_snprintf("%u", position);
		cv[cpIndex+4] = alleles[1];
		cv[cpIndex+5] = alleles[0];
		cv[cpIndex+6] = string_snprintf("%.8f", xfer.total / (2.0 * xfer.nrows));
	}
}

struct LoadDataPGENProvider2 : public LoadDataProvider2<LoadDataPGENProvider2> {
	int cpIndex;

	struct PgenFileInfoDtor {
		void operator()(PgenFileInfo *pfi) {
			PglErr err = kPglRetSuccess;
			CleanupPgfi(pfi, &err);
			if (err != kPglRetSuccess) mxThrow("CleanupPgfi not happy");
			if (pfi->vrtypes) aligned_free(pfi->vrtypes);
			delete pfi;
		}
	};
	typedef std::unique_ptr< PgenFileInfo, PgenFileInfoDtor > PgenFileInfoPtr;
	struct PgenReaderDtor {
		void operator()(PgenReader *pgr) {
			PglErr err = kPglRetSuccess;
			CleanupPgr(pgr, &err);
			if (err != kPglRetSuccess) mxThrow("CleanupPgr not happy");
			delete pgr;
		}
	};
	typedef std::unique_ptr< PgenReader, PgenReaderDtor > PgenReaderPtr;

	PgenFileInfoPtr pgen_info;
	PgenReaderPtr pgen_state;
  // unaligned ptr stored in [0], aligned stored in [1]
  PgrSampleSubsetIndex pssi;
  uintptr_t* pgen_subset_include_vec;
	uint32_t* pgen_subset_cumulative_popcounts;
	uintptr_t* pgen_genovec;
	uintptr_t* pgen_dosage_present;
	uint16_t* pgen_dosage_main;

	std::vector<int> srcDestMap;

	virtual const char *getName() { return "pgen"; };
	virtual void init(SEXP rObj);
	virtual void addCheckpointColumns(std::vector< std::string > &cp)
	{
		cpIndex = cp.size();
		cp.push_back("MAF");
	}
	virtual void loadRowImpl(int index);
	virtual int getNumVariants();
};

void LoadDataPGENProvider2::init(SEXP rObj)
{
	requireFile(rObj);

	srcDestMap.resize(srcRows);
	for (int rx=0, dr=0; rx < srcRows; ++rx) {
		if (skipRow(rx)) {
			srcDestMap[rx] = -1;
		} else {
			srcDestMap[rx] = dr++;
		}
	}
}

void myGenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct,
														std::function<bool(int)> skipFn, void* __restrict result)
{
  const uint64_t* table_alias = S_CAST(const uint64_t*, table16x8bx2);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
	int rx=0, dr=0;
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
					if (!skipFn(rx)) {
						memcpy(result_iter + dr, &(table_alias[(geno_word & 3) * 2]), 8);
					}
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
			if (!skipFn(rx)) {
				memcpy(result_iter + dr, &(table_alias[cur_2geno * 2]), 8);
				dr += 1;
			}
			rx += 1;
			if (!skipFn(rx)) {
				memcpy(result_iter + dr, &(table_alias[cur_2geno * 2])+1, 8);
				dr += 1;
			}
			rx += 1;
      geno_word >>= 4;
    }
  }
}

const double kGenoDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, NA_REAL);

void Dosage16ToDoubles(const uintptr_t* genoarr, const uintptr_t* dosage_present,
											 const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct,
											 std::function<bool(int)> skipFn,
											 std::vector<int> &srcDestMap, double* geno_double)
{
	const double* geno_double_pair_table = kGenoDoublePairs;
	myGenoarrLookup16x8bx2(genoarr, geno_double_pair_table, sample_ct, skipFn, geno_double);
  if (dosage_ct) {
		// Not all data has poor QC.
    const uint16_t* dosage_main_iter = dosage_main;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = dosage_present[0];
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
			int dr = srcDestMap[sample_uidx];
			if (dr >= 0) {
				geno_double[dr] = S_CAST(double, *dosage_main_iter) * 0.00006103515625;
			}
			dosage_main_iter += 1;
    }
  }
}

int LoadDataPGENProvider2::getNumVariants()
{ return pgen_info->raw_variant_ct; }

void LoadDataPGENProvider2::loadRowImpl(int index)
{
	if (columns.size() != 1) mxThrow("%s: pgen only has 1 column, not %d",
					 name, int(columns.size()));

	// adapted from plink-ng/2.0/Python/pgenlib.pyx
	if (!pgen_info) {
		pgen_info = PgenFileInfoPtr(new PgenFileInfo);
		PreinitPgfi(pgen_info.get());
		pgen_info->vrtypes = 0;
		uint32_t cur_variant_ct = 0xffffffffU;
		uint32_t cur_sample_ct = srcRows;
		PgenHeaderCtrl header_ctrl;
		uintptr_t pgfi_alloc_cacheline_ct;
		char errstr_buf[kPglErrstrBufBlen];
		if (PgfiInitPhase1(filePath.c_str(), cur_variant_ct, cur_sample_ct, 0, &header_ctrl,
				   pgen_info.get(), &pgfi_alloc_cacheline_ct, errstr_buf) != kPglRetSuccess) {
			mxThrow("%s: PgfiInitPhase1(%s) %s", name, filePath.c_str(), errstr_buf);
		}
		// No idea the purpose of these assertions, copied verbatim
		assert((header_ctrl & 0x30) == 0); // no alt allele counts
		assert((header_ctrl & 0xc0) != 0xc0); // no explicit nonref_flags
		if (pgen_info->raw_sample_ct == 0)
			mxThrow("%s: pgen file '%s' has no samples", name, filePath.c_str());
    unsigned char *pgfi_alloc; // is freed automatically by libpgen
		if (pgfi_alloc_cacheline_ct != 0) {
			if (cachealigned_malloc(pgfi_alloc_cacheline_ct * kCacheline, &pgfi_alloc))
				mxThrow("%s: cachealigned_malloc failed", name);
		}
		uint32_t max_vrec_width;
		uintptr_t pgr_alloc_cacheline_ct;
		if (PgfiInitPhase2(header_ctrl, 1, 1, 0, 0, pgen_info->raw_variant_ct,
				   &max_vrec_width, pgen_info.get(), pgfi_alloc, &pgr_alloc_cacheline_ct,
				   errstr_buf)) {
			mxThrow("%s: PgfiInitPhase2(%s) %s", name, filePath.c_str(), errstr_buf);
		}
		pgen_state = PgenReaderPtr(new PgenReader);
		PreinitPgr(pgen_state.get());
		uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * kCacheline;
		uint32_t file_sample_ct = pgen_info->raw_sample_ct;
		uintptr_t sample_subset_byte_ct = DivUp(file_sample_ct, kBitsPerVec) * kBytesPerVec;
		uintptr_t cumulative_popcounts_byte_ct =
			DivUp(file_sample_ct, kBitsPerWord * kInt32PerVec) * kBytesPerVec;
    uintptr_t genovec_byte_ct = DivUp(file_sample_ct, kNypsPerVec) * kBytesPerVec;
		uintptr_t dosage_main_byte_ct = DivUp(file_sample_ct, (2 * kInt32PerVec)) * kBytesPerVec;
		unsigned char* pgr_alloc;
		if (cachealigned_malloc(pgr_alloc_main_byte_ct +
                            (2 * kPglNypTransposeBatch + 5) * sample_subset_byte_ct +
                            cumulative_popcounts_byte_ct +
                            (1 + kPglNypTransposeBatch) * genovec_byte_ct +
                            dosage_main_byte_ct + kPglBitTransposeBufbytes +
                            4 * (kPglNypTransposeBatch * kPglNypTransposeBatch % 8), &pgr_alloc))
			mxThrow("%s: cachealigned_malloc failed", name);
		PglErr reterr = PgrInit(filePath.c_str(), max_vrec_width, pgen_info.get(),
					pgen_state.get(), pgr_alloc);
		if (reterr != kPglRetSuccess) {
			mxThrow("%s: PgrInit(%s) error code %d", name, filePath.c_str(), int(reterr));
		}

		unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
    pgen_subset_include_vec = (uintptr_t*)pgr_alloc_iter;
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		pgen_subset_cumulative_popcounts = (uint32_t*)pgr_alloc_iter;
		pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
		pgen_genovec = (uintptr_t*)pgr_alloc_iter;
		pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
		pgen_dosage_present = (uintptr_t*)pgr_alloc_iter;
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		pgen_dosage_main = (uint16_t*)pgr_alloc_iter;
		loadCounter += 1;
	}

	if (1+index > int(pgen_info->raw_variant_ct)) {
	  mxThrow("%s: out of data (record %d requested but only %d in file)",
		  name, 1+index, int(pgen_info->raw_variant_ct));
	}

	double maf = 0;
	if (colTypes[0] == COLUMNDATA_NUMERIC) {
		uint32_t dosage_ct;
		PglErr reterr = PgrGet1D(pgen_subset_include_vec, pssi,
					 pgen_info->raw_sample_ct, index, 1, pgen_state.get(), pgen_genovec,
					 pgen_dosage_present, pgen_dosage_main, &dosage_ct);
		if (reterr != kPglRetSuccess)
			mxThrow("%s: read_dosages(varient %d) error code %d", name, index, int(reterr));
		Dosage16ToDoubles(pgen_genovec, pgen_dosage_present, pgen_dosage_main,
											pgen_info->raw_sample_ct, dosage_ct, [&](int rx)->bool{ return skipRow(rx); },
											srcDestMap, stripeData[0].realData);

		double *rd = stripeData[0].realData;
		int nrows = 0;
		for (int rx=0; rx < destRows; ++rx) {
			if (!std::isfinite(rd[rx])) continue;
			maf += rd[rx];
			nrows += 1;
		}
		maf /= 2.0 * nrows;
	} else {
		stop("Treating genetic data as an ordinal factor is not implemented");
	}

	for (int cx=0; cx < int(columns.size()); ++cx) {
		(*rawCols)[ columns[cx] ].setBorrow(stripeData[cx]);
	}
	if (checkpoint) {
		auto &cv = *checkpointValues;
		cv[cpIndex] = string_snprintf("%.8f", maf);
	}
}

unsigned int DJBHash(const char *str, std::size_t len)
{
   unsigned int hash = 5381;

   for(std::size_t i = 0; i < len; i++) {
     hash = ((hash << 5) + hash) + str[i];
   }

   return hash;
}

void setup2(AddLoadDataProviderType aldp)
{
  std::size_t sz2[] = {
               sizeof(dataPtr),
               sizeof(LoadDataProviderBase2),
               sizeof(ColumnData)
  };
  auto apiHash = DJBHash((char*)sz2, sizeof(sz2));
	aldp(OPENMX_LOAD_DATA_API_VERSION, apiHash, new LoadDataPGENProvider2());
	aldp(OPENMX_LOAD_DATA_API_VERSION, apiHash, new LoadDataBGENProvider2());
}
