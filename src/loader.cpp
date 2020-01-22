#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include "pgenlib_internal.h"
#include "openmx.h"
#include "LoadDataAPI.h"

using namespace plink2;

struct BgenXfer {
	dataPtr dp;
	std::vector<double> prob;
	int row;
	BgenXfer(dataPtr &_dp) : dp(_dp) {};
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {}
	void set_min_max_ploidy(genfile::bgen::uint32_t min_ploidy, genfile::bgen::uint32_t max_ploidy,
				genfile::bgen::uint32_t min_entries, genfile::bgen::uint32_t max_entries)
	{
		row = 0;
		if (min_ploidy != 2 || max_ploidy != 2 || min_entries != 3 || max_entries != 3) {
			mxThrow("set_min_max_ploidy %u %u %u %u, not implemented",
				min_ploidy, max_ploidy, min_entries, max_entries);
		}
	}
	bool set_sample( std::size_t i ) { return true; }
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
			double dosage = prob[1] * 1 + prob[2] * 2;
			dp.realData[row++] = dosage;
		}
	}

	void set_value( genfile::bgen::uint32_t entry_i, genfile::MissingValue) {
		if (entry_i == 2) {
			dp.realData[row++] = NA_REAL;
		}
	}
};

class LoadDataBGENProvider : public LoadDataProvider<LoadDataBGENProvider> {
	int cpIndex;
	genfile::bgen::View::UniquePtr bgenView;

	virtual const char *getName() { return "bgen"; };
	virtual void init(SEXP rObj) {
		ProtectedSEXP Rbyrow(R_do_slot(rObj, Rf_install("byrow")));
		bool byrow = Rf_asLogical(Rbyrow);
		if (!byrow) mxThrow("byrow=FALSE is not implemented for bgen format");
		requireFile(rObj);
	}
	virtual void loadRowImpl(int index);
	virtual void addCheckpointColumns(std::vector< std::string > &cp)
	{
		cpIndex = cp.size();
		std::string c1 = "SNP";
		cp.push_back(c1);
		c1 = "RSID";
		cp.push_back(c1);
		c1 = "CHR";
		cp.push_back(c1);
		c1 = "BP";
		cp.push_back(c1);
		c1 = "A1";
		cp.push_back(c1);
		c1 = "A2";
		cp.push_back(c1);
	}
	virtual int getNumVariants();
};

int LoadDataBGENProvider::getNumVariants()
{
	if (bgenView.get() == 0) return 0;
	return bgenView->number_of_variants();
}

void LoadDataBGENProvider::loadRowImpl(int index)
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
		if (rows != int(bgenView->number_of_samples())) {
			mxThrow("%s: %s has %d rows but %s has %d samples",
				name, dataName, rows, filePath.c_str(),
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
	if (checkpoint) {
		auto &cv = *checkpointValues;
		cv[cpIndex] = SNPID;
		cv[cpIndex+1] = rsid;
		cv[cpIndex+2] = chromosome;
		cv[cpIndex+3] = string_snprintf("%u", position);
		cv[cpIndex+4] = alleles[1];
		cv[cpIndex+5] = alleles[0];
	}
	BgenXfer xfer(stripeData[0]);
	bgenView->read_genotype_data_block(xfer);
	curRecord += 1;

	auto &rc = *rawCols;
	for (int cx=0; cx < int(columns.size()); ++cx) {
		rc[ columns[cx] ].ptr = stripeData[cx];
	}
}

class LoadDataPGENProvider : public LoadDataProvider<LoadDataPGENProvider> {
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
	struct PgenReaderStructDtor {
		void operator()(PgenReaderStruct *pgr) {
			PglErr err = kPglRetSuccess;
			CleanupPgr(pgr, &err);
			if (err != kPglRetSuccess) mxThrow("CleanupPgr not happy");
			if (pgr->fread_buf) aligned_free(pgr->fread_buf);
			delete pgr;
		}
	};
	typedef std::unique_ptr< PgenReaderStruct, PgenReaderStructDtor > PgenReaderStructPtr;
		
	PgenFileInfoPtr pgen_info;
	PgenReaderStructPtr pgen_state;
	uintptr_t* pgen_subset_include_vec;
	uint32_t* pgen_subset_cumulative_popcounts;
	uintptr_t* pgen_genovec;
	uintptr_t* pgen_dosage_present;
	uint16_t* pgen_dosage_main;

	virtual const char *getName() { return "pgen"; };
	virtual void init(SEXP rObj) {
		requireFile(rObj);
	}
	virtual void loadRowImpl(int index);
	virtual int getNumVariants();
};

static const int kGenoToFactor[4] = {1, 2, 3, NA_INTEGER};

// TODO: investigate GenoarrLookup16x8bx2()
static void GenoarrToFactor(const uintptr_t* genoarr, uint32_t sample_ct, int *geno_out) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  int* write_iter = geno_out;
  uint32_t subgroup_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        return;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      *write_iter++ = kGenoToFactor[geno_word & 3];
      geno_word >>= 2;
    }
  }
}

void Dosage16ToDoubles(const double* geno_double_pair_table, const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double) {
  GenoarrLookup16x8bx2(genoarr, geno_double_pair_table, sample_ct, geno_double);
  if (dosage_ct) {
    const uint16_t* dosage_main_iter = dosage_main;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = dosage_present[0];
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
      // We use 2.0- to obtain exactly the same dosages as bgen format.
      // I guess it's arbitrary so why not?
      geno_double[sample_uidx] = 2.0 - S_CAST(double, *dosage_main_iter++) * 0.00006103515625;
    }
  }
}

const double kGenoDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, -9.0);

void Dosage16ToDoublesMinus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double)
{
	Dosage16ToDoubles(kGenoDoublePairs, genoarr, dosage_present, dosage_main, sample_ct, dosage_ct, geno_double);
}

int LoadDataPGENProvider::getNumVariants()
{ return pgen_info->raw_variant_ct; }

void LoadDataPGENProvider::loadRowImpl(int index)
{
	if (columns.size() != 1) mxThrow("%s: pgen only has 1 column, not %d",
					 name, int(columns.size()));

	// adapted from plink-ng/2.0/Python/pgenlib.pyx
	if (!pgen_info) {
		pgen_info = PgenFileInfoPtr(new PgenFileInfo);
		PreinitPgfi(pgen_info.get());
		pgen_info->vrtypes = 0;
		uint32_t cur_variant_ct = 0xffffffffU;
		uint32_t cur_sample_ct = rows;
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
		unsigned char* pgfi_alloc = 0;
		if (pgfi_alloc_cacheline_ct != 0) {
			if (cachealigned_malloc(pgfi_alloc_cacheline_ct * kCacheline, &pgfi_alloc))
				mxThrow("%s: cachealigned_malloc failed", name);
		}
		uint32_t max_vrec_width;
		uintptr_t pgr_alloc_cacheline_ct;
		if (PgfiInitPhase2(header_ctrl, 1, 1, 0, 0, pgen_info->raw_variant_ct,
				   &max_vrec_width, pgen_info.get(), pgfi_alloc, &pgr_alloc_cacheline_ct,
				   errstr_buf)) {
			if (pgfi_alloc && !pgen_info->vrtypes) aligned_free(pgfi_alloc);
			mxThrow("%s: PgfiInitPhase2(%s) %s", name, filePath.c_str(), errstr_buf);
		}
		pgen_state = PgenReaderStructPtr(new PgenReaderStruct);
		PreinitPgr(pgen_state.get());
		pgen_state->fread_buf = 0;
		uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * kCacheline;
		uint32_t file_sample_ct = pgen_info->raw_sample_ct;
		uintptr_t sample_subset_byte_ct = DivUp(file_sample_ct, kBitsPerVec) * kBytesPerVec;
		uintptr_t cumulative_popcounts_byte_ct =
			DivUp(file_sample_ct, kBitsPerWord * kInt32PerVec) * kBytesPerVec;
		uintptr_t genovec_byte_ct = DivUp(file_sample_ct, kQuatersPerVec) * kBytesPerVec;
		uintptr_t dosage_main_byte_ct = DivUp(file_sample_ct, (2 * kInt32PerVec)) * kBytesPerVec;
		unsigned char* pgr_alloc;
		if (cachealigned_malloc(pgr_alloc_main_byte_ct +
					(2 * kPglQuaterTransposeBatch + 5) * sample_subset_byte_ct +
					cumulative_popcounts_byte_ct +
					(1 + kPglQuaterTransposeBatch) * genovec_byte_ct +
					dosage_main_byte_ct, &pgr_alloc))
			mxThrow("%s: cachealigned_malloc failed", name);
		PglErr reterr = PgrInit(filePath.c_str(), max_vrec_width, pgen_info.get(),
					pgen_state.get(), pgr_alloc);
		if (reterr != kPglRetSuccess) {
			if (!pgen_state->fread_buf) aligned_free(pgr_alloc);
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

	if (colTypes[0] == COLUMNDATA_NUMERIC) {
		uint32_t dosage_ct;
		PglErr reterr = PgrGet1D(pgen_subset_include_vec, pgen_subset_cumulative_popcounts,
					 pgen_info->raw_sample_ct, index, 1, pgen_state.get(), pgen_genovec,
					 pgen_dosage_present, pgen_dosage_main, &dosage_ct);
		if (reterr != kPglRetSuccess)
			mxThrow("%s: read_dosages(varient %d) error code %d", name, index, int(reterr));
		Dosage16ToDoublesMinus9(pgen_genovec, pgen_dosage_present, pgen_dosage_main,
					pgen_info->raw_sample_ct, dosage_ct, stripeData[0].realData);
	} else {
		PglErr reterr = PgrGet1(pgen_subset_include_vec, pgen_subset_cumulative_popcounts,
					pgen_info->raw_sample_ct, index, 1, pgen_state.get(), pgen_genovec);
		if (reterr != kPglRetSuccess)
			mxThrow("%s: read(varient %d) error code %d", name, index, int(reterr));

		auto &rc = (*rawCols)[ columns[0] ];
		if (rc.levels.size() != 3) mxThrow("%s: pgen files contain data with 3 levels (not %d)",
						   name, int(rc.levels.size()));
		GenoarrToFactor(pgen_genovec, pgen_info->raw_sample_ct, stripeData[0].intData);
	}

	for (int cx=0; cx < int(columns.size()); ++cx) {
		(*rawCols)[ columns[cx] ].ptr = stripeData[cx];
	}
}

void setup2(AddLoadDataProviderType aldp)
{
	int sz = sizeof(LoadDataProviderBase);
	aldp(OPENMX_LOAD_DATA_API_VERSION, sz, new LoadDataPGENProvider());
	aldp(OPENMX_LOAD_DATA_API_VERSION, sz, new LoadDataBGENProvider());
}
