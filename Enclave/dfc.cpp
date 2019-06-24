/*********************************/
/* Author  - Byungkwon Choi      */
/* Contact - cbkbrad@kaist.ac.kr */
/*********************************/

/***********************************/
/* Modified by  - Juhyeng Han      */
/* Contact - sparkly9399@gmail.com */
/***********************************/

#include "dfc.h"

/*************************************************************************************/
#define INIT_HASH_SIZE 65536
#define RECURSIVE_BOUNDARY 5
/*************************************************************************************/

/*************************************************************************************/
#define PRINT_INFO
#define CT8_SWITCH

#define ENABLE_PROGRESSIVE_FILTERING
#define ENABLE_RECURSIVE

#define COUNT_MATCH
/*************************************************************************************/

/*************************************************************************************/
static unsigned char xlatcase[256];

/* For extracting position */
int pattern_interval;
int min_pattern_interval;

static int dfc_total_memory = 0;
static int dfc_pattern_memory = 0;
static int dfc_memory_dfs = (sizeof(uint8_t) * DF_SIZE_REAL) * 14;
static int dfc_memory_ct1 = 0;
static int dfc_memory_ct2 = 0;
static int dfc_memory_ct3 = 0;
static int dfc_memory_ct4 = 0;
static int dfc_memory_ct8 = 0;
/*************************************************************************************/

/*************************************************************************************/
static float my_sqrtf(float input, float x);
static void init_xlatcase();
static inline void ConvertCaseEx (unsigned char *d, unsigned char *s, int m);
static inline int my_strncmp(unsigned char *a, unsigned char *b, int n);
static inline int my_strncasecmp(unsigned char *a, unsigned char *b, int n);
static void * DFC_REALLOC(void *p, uint16_t n, dfcDataType type, dfcMemoryType type2 );
static void DFC_FREE(void *p, int n, dfcMemoryType type);
static void * DFC_MALLOC(int n, dfcMemoryType type );
static void Build_pattern(DFC_PATTERN *p, uint8_t *flag, uint8_t *temp, uint32_t i, int j, int k);
static inline DFC_PATTERN *DFC_InitHashLookup(DFC_STRUCTURE *ctx, uint8_t *pat, uint16_t patlen, int nocase);
static inline int DFC_InitHashAdd(DFC_STRUCTURE *ctx, DFC_PATTERN *p);
/*************************************************************************************/

static const uint64_t crc64_tab[256] = {
    UINT64_C(0x0000000000000000), UINT64_C(0x7ad870c830358979),
    UINT64_C(0xf5b0e190606b12f2), UINT64_C(0x8f689158505e9b8b),
    UINT64_C(0xc038e5739841b68f), UINT64_C(0xbae095bba8743ff6),
    UINT64_C(0x358804e3f82aa47d), UINT64_C(0x4f50742bc81f2d04),
    UINT64_C(0xab28ecb46814fe75), UINT64_C(0xd1f09c7c5821770c),
    UINT64_C(0x5e980d24087fec87), UINT64_C(0x24407dec384a65fe),
    UINT64_C(0x6b1009c7f05548fa), UINT64_C(0x11c8790fc060c183),
    UINT64_C(0x9ea0e857903e5a08), UINT64_C(0xe478989fa00bd371),
    UINT64_C(0x7d08ff3b88be6f81), UINT64_C(0x07d08ff3b88be6f8),
    UINT64_C(0x88b81eabe8d57d73), UINT64_C(0xf2606e63d8e0f40a),
    UINT64_C(0xbd301a4810ffd90e), UINT64_C(0xc7e86a8020ca5077),
    UINT64_C(0x4880fbd87094cbfc), UINT64_C(0x32588b1040a14285),
    UINT64_C(0xd620138fe0aa91f4), UINT64_C(0xacf86347d09f188d),
    UINT64_C(0x2390f21f80c18306), UINT64_C(0x594882d7b0f40a7f),
    UINT64_C(0x1618f6fc78eb277b), UINT64_C(0x6cc0863448deae02),
    UINT64_C(0xe3a8176c18803589), UINT64_C(0x997067a428b5bcf0),
    UINT64_C(0xfa11fe77117cdf02), UINT64_C(0x80c98ebf2149567b),
    UINT64_C(0x0fa11fe77117cdf0), UINT64_C(0x75796f2f41224489),
    UINT64_C(0x3a291b04893d698d), UINT64_C(0x40f16bccb908e0f4),
    UINT64_C(0xcf99fa94e9567b7f), UINT64_C(0xb5418a5cd963f206),
    UINT64_C(0x513912c379682177), UINT64_C(0x2be1620b495da80e),
    UINT64_C(0xa489f35319033385), UINT64_C(0xde51839b2936bafc),
    UINT64_C(0x9101f7b0e12997f8), UINT64_C(0xebd98778d11c1e81),
    UINT64_C(0x64b116208142850a), UINT64_C(0x1e6966e8b1770c73),
    UINT64_C(0x8719014c99c2b083), UINT64_C(0xfdc17184a9f739fa),
    UINT64_C(0x72a9e0dcf9a9a271), UINT64_C(0x08719014c99c2b08),
    UINT64_C(0x4721e43f0183060c), UINT64_C(0x3df994f731b68f75),
    UINT64_C(0xb29105af61e814fe), UINT64_C(0xc849756751dd9d87),
    UINT64_C(0x2c31edf8f1d64ef6), UINT64_C(0x56e99d30c1e3c78f),
    UINT64_C(0xd9810c6891bd5c04), UINT64_C(0xa3597ca0a188d57d),
    UINT64_C(0xec09088b6997f879), UINT64_C(0x96d1784359a27100),
    UINT64_C(0x19b9e91b09fcea8b), UINT64_C(0x636199d339c963f2),
    UINT64_C(0xdf7adabd7a6e2d6f), UINT64_C(0xa5a2aa754a5ba416),
    UINT64_C(0x2aca3b2d1a053f9d), UINT64_C(0x50124be52a30b6e4),
    UINT64_C(0x1f423fcee22f9be0), UINT64_C(0x659a4f06d21a1299),
    UINT64_C(0xeaf2de5e82448912), UINT64_C(0x902aae96b271006b),
    UINT64_C(0x74523609127ad31a), UINT64_C(0x0e8a46c1224f5a63),
    UINT64_C(0x81e2d7997211c1e8), UINT64_C(0xfb3aa75142244891),
    UINT64_C(0xb46ad37a8a3b6595), UINT64_C(0xceb2a3b2ba0eecec),
    UINT64_C(0x41da32eaea507767), UINT64_C(0x3b024222da65fe1e),
    UINT64_C(0xa2722586f2d042ee), UINT64_C(0xd8aa554ec2e5cb97),
    UINT64_C(0x57c2c41692bb501c), UINT64_C(0x2d1ab4dea28ed965),
    UINT64_C(0x624ac0f56a91f461), UINT64_C(0x1892b03d5aa47d18),
    UINT64_C(0x97fa21650afae693), UINT64_C(0xed2251ad3acf6fea),
    UINT64_C(0x095ac9329ac4bc9b), UINT64_C(0x7382b9faaaf135e2),
    UINT64_C(0xfcea28a2faafae69), UINT64_C(0x8632586aca9a2710),
    UINT64_C(0xc9622c4102850a14), UINT64_C(0xb3ba5c8932b0836d),
    UINT64_C(0x3cd2cdd162ee18e6), UINT64_C(0x460abd1952db919f),
    UINT64_C(0x256b24ca6b12f26d), UINT64_C(0x5fb354025b277b14),
    UINT64_C(0xd0dbc55a0b79e09f), UINT64_C(0xaa03b5923b4c69e6),
    UINT64_C(0xe553c1b9f35344e2), UINT64_C(0x9f8bb171c366cd9b),
    UINT64_C(0x10e3202993385610), UINT64_C(0x6a3b50e1a30ddf69),
    UINT64_C(0x8e43c87e03060c18), UINT64_C(0xf49bb8b633338561),
    UINT64_C(0x7bf329ee636d1eea), UINT64_C(0x012b592653589793),
    UINT64_C(0x4e7b2d0d9b47ba97), UINT64_C(0x34a35dc5ab7233ee),
    UINT64_C(0xbbcbcc9dfb2ca865), UINT64_C(0xc113bc55cb19211c),
    UINT64_C(0x5863dbf1e3ac9dec), UINT64_C(0x22bbab39d3991495),
    UINT64_C(0xadd33a6183c78f1e), UINT64_C(0xd70b4aa9b3f20667),
    UINT64_C(0x985b3e827bed2b63), UINT64_C(0xe2834e4a4bd8a21a),
    UINT64_C(0x6debdf121b863991), UINT64_C(0x1733afda2bb3b0e8),
    UINT64_C(0xf34b37458bb86399), UINT64_C(0x8993478dbb8deae0),
    UINT64_C(0x06fbd6d5ebd3716b), UINT64_C(0x7c23a61ddbe6f812),
    UINT64_C(0x3373d23613f9d516), UINT64_C(0x49aba2fe23cc5c6f),
    UINT64_C(0xc6c333a67392c7e4), UINT64_C(0xbc1b436e43a74e9d),
    UINT64_C(0x95ac9329ac4bc9b5), UINT64_C(0xef74e3e19c7e40cc),
    UINT64_C(0x601c72b9cc20db47), UINT64_C(0x1ac40271fc15523e),
    UINT64_C(0x5594765a340a7f3a), UINT64_C(0x2f4c0692043ff643),
    UINT64_C(0xa02497ca54616dc8), UINT64_C(0xdafce7026454e4b1),
    UINT64_C(0x3e847f9dc45f37c0), UINT64_C(0x445c0f55f46abeb9),
    UINT64_C(0xcb349e0da4342532), UINT64_C(0xb1eceec59401ac4b),
    UINT64_C(0xfebc9aee5c1e814f), UINT64_C(0x8464ea266c2b0836),
    UINT64_C(0x0b0c7b7e3c7593bd), UINT64_C(0x71d40bb60c401ac4),
    UINT64_C(0xe8a46c1224f5a634), UINT64_C(0x927c1cda14c02f4d),
    UINT64_C(0x1d148d82449eb4c6), UINT64_C(0x67ccfd4a74ab3dbf),
    UINT64_C(0x289c8961bcb410bb), UINT64_C(0x5244f9a98c8199c2),
    UINT64_C(0xdd2c68f1dcdf0249), UINT64_C(0xa7f41839ecea8b30),
    UINT64_C(0x438c80a64ce15841), UINT64_C(0x3954f06e7cd4d138),
    UINT64_C(0xb63c61362c8a4ab3), UINT64_C(0xcce411fe1cbfc3ca),
    UINT64_C(0x83b465d5d4a0eece), UINT64_C(0xf96c151de49567b7),
    UINT64_C(0x76048445b4cbfc3c), UINT64_C(0x0cdcf48d84fe7545),
    UINT64_C(0x6fbd6d5ebd3716b7), UINT64_C(0x15651d968d029fce),
    UINT64_C(0x9a0d8ccedd5c0445), UINT64_C(0xe0d5fc06ed698d3c),
    UINT64_C(0xaf85882d2576a038), UINT64_C(0xd55df8e515432941),
    UINT64_C(0x5a3569bd451db2ca), UINT64_C(0x20ed197575283bb3),
    UINT64_C(0xc49581ead523e8c2), UINT64_C(0xbe4df122e51661bb),
    UINT64_C(0x3125607ab548fa30), UINT64_C(0x4bfd10b2857d7349),
    UINT64_C(0x04ad64994d625e4d), UINT64_C(0x7e7514517d57d734),
    UINT64_C(0xf11d85092d094cbf), UINT64_C(0x8bc5f5c11d3cc5c6),
    UINT64_C(0x12b5926535897936), UINT64_C(0x686de2ad05bcf04f),
    UINT64_C(0xe70573f555e26bc4), UINT64_C(0x9ddd033d65d7e2bd),
    UINT64_C(0xd28d7716adc8cfb9), UINT64_C(0xa85507de9dfd46c0),
    UINT64_C(0x273d9686cda3dd4b), UINT64_C(0x5de5e64efd965432),
    UINT64_C(0xb99d7ed15d9d8743), UINT64_C(0xc3450e196da80e3a),
    UINT64_C(0x4c2d9f413df695b1), UINT64_C(0x36f5ef890dc31cc8),
    UINT64_C(0x79a59ba2c5dc31cc), UINT64_C(0x037deb6af5e9b8b5),
    UINT64_C(0x8c157a32a5b7233e), UINT64_C(0xf6cd0afa9582aa47),
    UINT64_C(0x4ad64994d625e4da), UINT64_C(0x300e395ce6106da3),
    UINT64_C(0xbf66a804b64ef628), UINT64_C(0xc5bed8cc867b7f51),
    UINT64_C(0x8aeeace74e645255), UINT64_C(0xf036dc2f7e51db2c),
    UINT64_C(0x7f5e4d772e0f40a7), UINT64_C(0x05863dbf1e3ac9de),
    UINT64_C(0xe1fea520be311aaf), UINT64_C(0x9b26d5e88e0493d6),
    UINT64_C(0x144e44b0de5a085d), UINT64_C(0x6e963478ee6f8124),
    UINT64_C(0x21c640532670ac20), UINT64_C(0x5b1e309b16452559),
    UINT64_C(0xd476a1c3461bbed2), UINT64_C(0xaeaed10b762e37ab),
    UINT64_C(0x37deb6af5e9b8b5b), UINT64_C(0x4d06c6676eae0222),
    UINT64_C(0xc26e573f3ef099a9), UINT64_C(0xb8b627f70ec510d0),
    UINT64_C(0xf7e653dcc6da3dd4), UINT64_C(0x8d3e2314f6efb4ad),
    UINT64_C(0x0256b24ca6b12f26), UINT64_C(0x788ec2849684a65f),
    UINT64_C(0x9cf65a1b368f752e), UINT64_C(0xe62e2ad306bafc57),
    UINT64_C(0x6946bb8b56e467dc), UINT64_C(0x139ecb4366d1eea5),
    UINT64_C(0x5ccebf68aecec3a1), UINT64_C(0x2616cfa09efb4ad8),
    UINT64_C(0xa97e5ef8cea5d153), UINT64_C(0xd3a62e30fe90582a),
    UINT64_C(0xb0c7b7e3c7593bd8), UINT64_C(0xca1fc72bf76cb2a1),
    UINT64_C(0x45775673a732292a), UINT64_C(0x3faf26bb9707a053),
    UINT64_C(0x70ff52905f188d57), UINT64_C(0x0a2722586f2d042e),
    UINT64_C(0x854fb3003f739fa5), UINT64_C(0xff97c3c80f4616dc),
    UINT64_C(0x1bef5b57af4dc5ad), UINT64_C(0x61372b9f9f784cd4),
    UINT64_C(0xee5fbac7cf26d75f), UINT64_C(0x9487ca0fff135e26),
    UINT64_C(0xdbd7be24370c7322), UINT64_C(0xa10fceec0739fa5b),
    UINT64_C(0x2e675fb4576761d0), UINT64_C(0x54bf2f7c6752e8a9),
    UINT64_C(0xcdcf48d84fe75459), UINT64_C(0xb71738107fd2dd20),
    UINT64_C(0x387fa9482f8c46ab), UINT64_C(0x42a7d9801fb9cfd2),
    UINT64_C(0x0df7adabd7a6e2d6), UINT64_C(0x772fdd63e7936baf),
    UINT64_C(0xf8474c3bb7cdf024), UINT64_C(0x829f3cf387f8795d),
    UINT64_C(0x66e7a46c27f3aa2c), UINT64_C(0x1c3fd4a417c62355),
    UINT64_C(0x935745fc4798b8de), UINT64_C(0xe98f353477ad31a7),
    UINT64_C(0xa6df411fbfb21ca3), UINT64_C(0xdc0731d78f8795da),
    UINT64_C(0x536fa08fdfd90e51), UINT64_C(0x29b7d047efec8728),
};

uint64_t crc64(uint64_t crc, uint64_t s, uint64_t l)
{
    uint64_t j;
    uint64_t mask = 0x00000000000000FF;
    for (j = 0; j < l; j++) {
        uint8_t byte = (s & mask) >> (j * 8);
        mask = mask << (j * 8);
        crc = crc64_tab[(uint8_t)crc ^ byte] ^ (crc >> 8);
    }
    return crc;
}

/*
*  Create a new DFC state machine
*/
DFC_STRUCTURE * DFC_New (void)
{
    DFC_STRUCTURE * p;

    init_xlatcase ();

    dfc_total_memory = 0;
    dfc_pattern_memory = 0;

    p = (DFC_STRUCTURE *)DFC_MALLOC(sizeof (DFC_STRUCTURE), DFC_MEMORY_TYPE__DFC);
    MEMASSERT_DFC (p, "DFC_New");

    if (p)
    {
        memset (p, 0, sizeof (DFC_STRUCTURE));
        p->init_hash = (DFC_PATTERN**)malloc(sizeof(DFC_PATTERN *) * INIT_HASH_SIZE);
        if(p->init_hash == NULL){
            sgx_exit(1);
        }
        memset(p->init_hash, 0, sizeof(DFC_PATTERN *) * INIT_HASH_SIZE);
    }

    return p;
}

static void DFC_FreePattern(DFC_PATTERN *p){
    if(p->patrn != NULL){
        free(p->patrn);
    }

    if(p->casepatrn != NULL){
        free(p->casepatrn);
    }

    return;
}

static void DFC_FreePatternList(DFC_STRUCTURE *dfc){
    DFC_PATTERN *plist;
    DFC_PATTERN *p_next;

	for(plist = dfc->dfcPatterns; plist != NULL;){
        DFC_FreePattern(plist);
        p_next = plist->next;
        free(plist);
        plist = p_next;
    }

    return;
}

void DFC_FreeStructure(DFC_STRUCTURE *dfc){
    BUC_CNT_TYPE j, l;
    int i,k;
    if(dfc == NULL)
        return;

    if(dfc->dfcPatterns != NULL){
        DFC_FreePatternList(dfc);
    }

    if(dfc->dfcMatchList != NULL){
        free(dfc->dfcMatchList);
    }

    for(i = 0; i < CT2_TABLE_SIZE; i++){
        for(j = 0; j < dfc->CompactTable2[i].cnt; j++){
            free(dfc->CompactTable2[i].array[j].pid);

            if(dfc->CompactTable2[i].array[j].DirectFilter != NULL){
                free(dfc->CompactTable2[i].array[j].DirectFilter);
            }

            if(dfc->CompactTable2[i].array[j].CompactTable != NULL){
                for(k = 0; k < RECURSIVE_CT_SIZE; k++){
                    for(l = 0; l < dfc->CompactTable2[i].array[j].CompactTable[k].cnt; l++){
                        free(dfc->CompactTable2[i].array[j].CompactTable[k].array[l].pid);
                    }
                    free(dfc->CompactTable2[i].array[j].CompactTable[k].array);
                }
                free(dfc->CompactTable2[i].array[j].CompactTable);
            }
        
        }
    }

    for(i = 0; i < CT4_TABLE_SIZE; i++){
        for(j = 0; j < dfc->CompactTable4[i].cnt; j++){
            free(dfc->CompactTable4[i].array[j].pid);

            if(dfc->CompactTable4[i].array[j].DirectFilter != NULL){
                free(dfc->CompactTable4[i].array[j].DirectFilter);
            }

            if(dfc->CompactTable4[i].array[j].CompactTable != NULL){
                for(k = 0; k < RECURSIVE_CT_SIZE; k++){
                    for(l = 0; l < dfc->CompactTable4[i].array[j].CompactTable[k].cnt; l++){
                        free(dfc->CompactTable4[i].array[j].CompactTable[k].array[l].pid);
                    }
                    free(dfc->CompactTable4[i].array[j].CompactTable[k].array);
                }
                free(dfc->CompactTable4[i].array[j].CompactTable);
            }
        
        }
    }

    for(i = 0; i < CT8_TABLE_SIZE; i++){
        for(j = 0; j < dfc->CompactTable8[i].cnt; j++){
            free(dfc->CompactTable8[i].array[j].pid);

            if(dfc->CompactTable8[i].array[j].DirectFilter != NULL){
                free(dfc->CompactTable8[i].array[j].DirectFilter);
            }

            if(dfc->CompactTable8[i].array[j].CompactTable != NULL){
                for(k = 0; k < RECURSIVE_CT_SIZE; k++){
                    for(l = 0; l < dfc->CompactTable8[i].array[j].CompactTable[k].cnt; l++){
                        free(dfc->CompactTable8[i].array[j].CompactTable[k].array[l].pid);
                    }
                    free(dfc->CompactTable8[i].array[j].CompactTable[k].array);
                }
                free(dfc->CompactTable8[i].array[j].CompactTable);
            }
        
        }
    }

    free(dfc);
}


/*
*  Add a pattern to the list of patterns
*
*
* \param dfc    Pointer to the DFC structure
* \param pat    Pointer to the pattern
* \param n      Pattern length
* \param nocase Flag for case-sensitivity (0 means case-insensitive, 1 means the opposite)
* \param sid    External id
*
* \retval   0 On success to add new pattern.
* \retval   1 On success to add sid.
*/
int DFC_AddPattern (DFC_STRUCTURE * dfc, unsigned char *pat, int n, int nocase, PID_TYPE sid)
{
    DFC_PATTERN * plist = DFC_InitHashLookup(dfc, pat, n, nocase);

    if(plist == NULL){
        plist = (DFC_PATTERN *) DFC_MALLOC(sizeof (DFC_PATTERN), DFC_MEMORY_TYPE__PATTERN);
        MEMASSERT_DFC (plist, "DFC_AddPattern");
        memset(plist, 0, sizeof(DFC_PATTERN));

        plist->patrn = (unsigned char *)DFC_MALLOC(n, DFC_MEMORY_TYPE__PATTERN);
        MEMASSERT_DFC (plist->patrn, "DFC_AddPattern");

        ConvertCaseEx(plist->patrn, pat, n);

        plist->casepatrn = (unsigned char *)DFC_MALLOC(n,DFC_MEMORY_TYPE__PATTERN);
        MEMASSERT_DFC (plist->casepatrn, "DFC_AddPattern");

        memcpy (plist->casepatrn, pat, n);

        plist->n      = n;
        plist->nocase = nocase;
        plist->iid    = dfc->numPatterns; // internal id
        plist->next   = NULL;

        DFC_InitHashAdd(dfc, plist);

        /* sid update */
        plist->sids_size = 1;
        plist->sids = (PID_TYPE *) DFC_MALLOC(sizeof(PID_TYPE), DFC_MEMORY_TYPE__PATTERN);
        MEMASSERT_DFC (plist->sids, "DFC_AddPattern");
        plist->sids[0] = sid;

        /* Add this pattern to the list */
        dfc->numPatterns++;

        return 0;
    }else{
        int found = 0;
        uint32_t x = 0;

        for (x = 0; x < plist->sids_size; x++) {
            if (plist->sids[x] == sid) {
                found = 1;
                break;
            }
        }

        if (!found) {
            PID_TYPE *sids = (PID_TYPE *)DFC_REALLOC(plist->sids, plist->sids_size + 1,
                                                     DFC_PID_TYPE, DFC_MEMORY_TYPE__PATTERN);
            plist->sids = sids;
            plist->sids[plist->sids_size] = sid;
            plist->sids_size++;
        }

        return 1;
    }
}


void DFC_PrintInfo(DFC_STRUCTURE* dfc){
	DFC_PATTERN *plist;
    BUC_CNT_TYPE j;
	int i;
	int nl1 = 0, nl2 = 0, nl3 = 0, nl4 = 0, nl8 = 0;

	int nocase_pat_cnt = 0;
	for(plist = dfc->dfcPatterns; plist != NULL; plist = plist->next){
        if (plist->n == 1) {
			nl1 ++;
        }else if (plist->n == 2) {
			nl2 ++;
        }else if (plist->n == 3) {
			nl3 ++;
        }else if (plist->n >= 4 && plist->n < 8) {
			nl4 ++;
        }else if (plist->n >= 8) {
			nl8 ++;
		}

		if(plist->nocase){
			nocase_pat_cnt++;
		}
	}

	int nb1 = 0, nb2 = 0,nb3 = 0,nb4 = 0,nb5 = 0, nb8 = 0;
	int nb9 = 0;
	for(i = 0; i < DF_SIZE_REAL ;i++)
		nb1 += __builtin_popcount(dfc->DirectFilter1[i]);
	for(i = 0; i < DF_SIZE_REAL ;i++)
		nb2 += __builtin_popcount(dfc->ADD_DF_4_plus[i]);
	for(i = 0; i < DF_SIZE_REAL ;i++)
		nb3 += __builtin_popcount(dfc->ADD_DF_8_1[i]);
	for(i = 0; i < 256 ;i++)
		nb4 += __builtin_popcount(dfc->cDF0[i]);
	for(i = 0; i < DF_SIZE_REAL ;i++)
		nb8 += __builtin_popcount(dfc->ADD_DF_4_1[i]);
	for(i = 0; i < DF_SIZE_REAL ;i++)
		nb9 += __builtin_popcount(dfc->cDF2[i]);

	/* CT 1*/
	int ct1_array_cnt = 0;
	int ct1_pid_cnt = 0;
	float ct1_pid_std_dev = 0;
	for(i=0; i < CT1_TABLE_SIZE; i++){
		ct1_pid_cnt += dfc->CompactTable1[i].cnt;
		if(dfc->CompactTable1[i].cnt != 0)
			ct1_array_cnt++;
	}
	for(i=0; i < CT1_TABLE_SIZE; i++){
		ct1_pid_std_dev += (dfc->CompactTable1[i].cnt - ((float)ct1_pid_cnt/CT1_TABLE_SIZE)) *
								(dfc->CompactTable1[i].cnt - ((float)ct1_pid_cnt/CT1_TABLE_SIZE));
	}
	ct1_pid_std_dev /= ct1_pid_cnt;
	ct1_pid_std_dev = my_sqrtf(ct1_pid_std_dev,ct1_pid_std_dev);

	/* CT 2*/
	BUC_CNT_TYPE ct2_array_min = 9999999;
	BUC_CNT_TYPE ct2_array_max = 0;
	PID_CNT_TYPE ct2_pid_min = 9999999;
	PID_CNT_TYPE ct2_pid_max = 0;
	int ct2_array_cnt = 0;
	int ct2_array_tot_cnt = 0;
	int ct2_pid_cnt = 0;
	int ct2_pid_tot_cnt = 0;
	float ct2_array_std_dev = 0;
	float ct2_pid_std_dev = 0;
	for(i=0; i < CT2_TABLE_SIZE; i++){
		ct2_array_tot_cnt += dfc->CompactTable2[i].cnt;	
		if(ct2_array_min > dfc->CompactTable2[i].cnt) ct2_array_min = dfc->CompactTable2[i].cnt;
		if(ct2_array_max < dfc->CompactTable2[i].cnt) ct2_array_max = dfc->CompactTable2[i].cnt;
		if(dfc->CompactTable2[i].cnt != 0)
			ct2_array_cnt++;
		for(j = 0; j < dfc->CompactTable2[i].cnt; j++){
			ct2_pid_tot_cnt += dfc->CompactTable2[i].array[j].cnt;
			if(ct2_pid_min > dfc->CompactTable2[i].array[j].cnt) 
                ct2_pid_min = dfc->CompactTable2[i].array[j].cnt;
			if(ct2_pid_max < dfc->CompactTable2[i].array[j].cnt) 
                ct2_pid_max = dfc->CompactTable2[i].array[j].cnt;
			if(dfc->CompactTable2[i].array[j].cnt != 0)
				ct2_pid_cnt++;
		}
	}
	for(i=0; i < CT2_TABLE_SIZE; i++){
		if(dfc->CompactTable2[i].cnt != 0){
			ct2_array_std_dev += (dfc->CompactTable2[i].cnt - ((float)ct2_array_tot_cnt/ct2_array_cnt))*
								 (dfc->CompactTable2[i].cnt - ((float)ct2_array_tot_cnt/ct2_array_cnt));	
		}
		for(j = 0; j < dfc->CompactTable2[i].cnt; j++){
			if(dfc->CompactTable2[i].array[j].cnt != 0){
				ct2_pid_std_dev += (dfc->CompactTable2[i].array[j].cnt 
                                    - ((float)ct2_pid_tot_cnt/ct2_pid_cnt)) *
								   (dfc->CompactTable2[i].array[j].cnt 
                                    - ((float)ct2_pid_tot_cnt/ct2_pid_cnt));
			}
		}
	}

	ct2_array_std_dev /=ct2_array_cnt;
	ct2_pid_std_dev /=ct2_pid_cnt;

	ct2_array_std_dev = my_sqrtf(ct2_array_std_dev, ct2_array_std_dev);
	ct2_pid_std_dev = my_sqrtf(ct2_pid_std_dev, ct2_pid_std_dev);

	/* CT4 */
	BUC_CNT_TYPE ct4_array_min = 9999999;
	BUC_CNT_TYPE ct4_array_max = 0;
	PID_CNT_TYPE ct4_pid_min = 9999999;
	PID_CNT_TYPE ct4_pid_max = 0;
	int ct4_array_cnt = 0;
	int ct4_mt2f_array_cnt = 0; // # of buckets that has more than 2 frags (mt2f : more than 2 fragments)
	int ct4_mt2f_frag_cnt = 0; 
	int ct4_array_tot_cnt = 0;
	int ct4_pid_cnt = 0;
	int ct4_pid_tot_cnt = 0;
	float ct4_array_std_dev = 0;
	float ct4_pid_std_dev = 0;
	for(i=0; i < CT4_TABLE_SIZE; i++){
		ct4_array_tot_cnt += dfc->CompactTable4[i].cnt;	
		if(ct4_array_min > dfc->CompactTable4[i].cnt) ct4_array_min = dfc->CompactTable4[i].cnt;
		if(ct4_array_max < dfc->CompactTable4[i].cnt) ct4_array_max = dfc->CompactTable4[i].cnt;
		if(dfc->CompactTable4[i].cnt != 0)
			ct4_array_cnt++;
		if(dfc->CompactTable4[i].cnt >= 2){
			ct4_mt2f_array_cnt++;
            ct4_mt2f_frag_cnt += (dfc->CompactTable4[i].cnt - 1);
        }
		for(j = 0; j < dfc->CompactTable4[i].cnt; j++){
			ct4_pid_tot_cnt += dfc->CompactTable4[i].array[j].cnt;
			if(ct4_pid_min > dfc->CompactTable4[i].array[j].cnt) 
                ct4_pid_min = dfc->CompactTable4[i].array[j].cnt;
			if(ct4_pid_max < dfc->CompactTable4[i].array[j].cnt) 
                ct4_pid_max = dfc->CompactTable4[i].array[j].cnt;
			if(dfc->CompactTable4[i].array[j].cnt != 0)
				ct4_pid_cnt++;
		}
	}
	for(i=0; i < CT4_TABLE_SIZE; i++){
		if(dfc->CompactTable4[i].cnt != 0){
			ct4_array_std_dev += (dfc->CompactTable4[i].cnt - ((float)ct4_array_tot_cnt/ct4_array_cnt))*
								 (dfc->CompactTable4[i].cnt - ((float)ct4_array_tot_cnt/ct4_array_cnt));	
		}
		for(j = 0; j < dfc->CompactTable4[i].cnt; j++){
			if(dfc->CompactTable4[i].array[j].cnt != 0){
				ct4_pid_std_dev += (dfc->CompactTable4[i].array[j].cnt 
                                    - ((float)ct4_pid_tot_cnt/ct4_pid_cnt)) *
								   (dfc->CompactTable4[i].array[j].cnt 
                                    - ((float)ct4_pid_tot_cnt/ct4_pid_cnt));
			}
		}
	}

	ct4_array_std_dev /=ct4_array_cnt;
	ct4_pid_std_dev /=ct4_pid_cnt;

	ct4_array_std_dev = my_sqrtf(ct4_array_std_dev, ct4_array_std_dev);
	ct4_pid_std_dev = my_sqrtf(ct4_pid_std_dev, ct4_pid_std_dev);

	/* CT8 */
	BUC_CNT_TYPE ct8_array_min = 9999999;
	BUC_CNT_TYPE ct8_array_max = 0;
	PID_CNT_TYPE ct8_pid_min = 9999999;
	PID_CNT_TYPE ct8_pid_max = 0;
	int ct8_array_cnt = 0;
    int ct8_mt2f_array_cnt = 0; // # of buckets that has more than 2 frags (mt2f : more than 2 fragments)
    int ct8_mt2f_frag_cnt = 0;
	int ct8_array_tot_cnt = 0;
	int ct8_pid_cnt = 0;
	int ct8_pid_tot_cnt = 0;
	float ct8_array_std_dev = 0;
	float ct8_pid_std_dev = 0;
	for(i=0; i < CT8_TABLE_SIZE; i++){
		ct8_array_tot_cnt += dfc->CompactTable8[i].cnt;	
		if(ct8_array_min > dfc->CompactTable8[i].cnt) ct8_array_min = dfc->CompactTable8[i].cnt;
		if(ct8_array_max < dfc->CompactTable8[i].cnt) ct8_array_max = dfc->CompactTable8[i].cnt;
		if(dfc->CompactTable8[i].cnt != 0)
			ct8_array_cnt++;
		if(dfc->CompactTable8[i].cnt >= 2){
			ct8_mt2f_array_cnt++;
            ct8_mt2f_frag_cnt += (dfc->CompactTable8[i].cnt - 1);
        }
		for(j = 0; j < dfc->CompactTable8[i].cnt; j++){
			ct8_pid_tot_cnt += dfc->CompactTable8[i].array[j].cnt;
			if(ct8_pid_min > dfc->CompactTable8[i].array[j].cnt) 
                ct8_pid_min = dfc->CompactTable8[i].array[j].cnt;
			if(ct8_pid_max < dfc->CompactTable8[i].array[j].cnt) 
                ct8_pid_max = dfc->CompactTable8[i].array[j].cnt;
			if(dfc->CompactTable8[i].array[j].cnt != 0)
				ct8_pid_cnt++;
		}
	}
	for(i=0; i < CT8_TABLE_SIZE; i++){
		if(dfc->CompactTable8[i].cnt != 0){
			ct8_array_std_dev += (dfc->CompactTable8[i].cnt - ((float)ct8_array_tot_cnt/ct8_array_cnt))*
								 (dfc->CompactTable8[i].cnt - ((float)ct8_array_tot_cnt/ct8_array_cnt));	
		}
		for(j = 0; j < dfc->CompactTable8[i].cnt; j++){
			if(dfc->CompactTable8[i].array[j].cnt != 0){
				ct8_pid_std_dev += (dfc->CompactTable8[i].array[j].cnt 
                                    - ((float)ct8_pid_tot_cnt/ct8_pid_cnt)) *
								   (dfc->CompactTable8[i].array[j].cnt 
                                    - ((float)ct8_pid_tot_cnt/ct8_pid_cnt));
			}
		}
	}

	ct8_array_std_dev /=ct8_array_cnt;
	ct8_pid_std_dev /=ct8_pid_cnt;

	ct8_array_std_dev = my_sqrtf(ct8_array_std_dev, ct8_array_std_dev);
	ct8_pid_std_dev = my_sqrtf(ct8_pid_std_dev, ct8_pid_std_dev);

	printf("\n");
	printf("+- [ Direct Filter + Compact Table(DFC) Summary ] -------------------------------------\n");
	printf("| Patterns: %5d (Case-Insensitive patterns: %d)\n", dfc->numPatterns, nocase_pat_cnt);
	printf("|   - 1B      : %5d\n", nl1);
	printf("|   - 2B      : %5d\n", nl2);
	printf("|   - 3B      : %5d\n", nl3);
	printf("|   - 4B ~ 7B : %5d\n", nl4);
	printf("|   - 8B ~    : %5d\n", nl8);
	printf("|\n");
	printf("| < Direct Filter Density > \n");
	printf("|   1. All patterns(DF1)            : %.6f (%4d)\n", (double)nb1/DF_SIZE, nb1);
	printf("|     1) 1B patterns(DF4)           : %.6f (%4d)\n", (double)nb4/DF_SIZE, nb4);
	printf("|     2) 2B patterns(DF5)           : %.6f (%4d)\n", (double)nb5/DF_SIZE, nb5);
	printf("|     4) 4B ~ patterns(DF2)         : %.6f (%4d)\n", (double)nb2/DF_SIZE, nb2);
	printf("|       (1) 4B ~ 7B patterns(DF8)   : %.6f (%4d)\n", (double)nb8/DF_SIZE, nb8);
	printf("|          I. 4B ~ 7B patterns(DF9) : %.6f (%4d)\n", (double)nb9/DF_SIZE, nb9);
	printf("|       (2) 8B ~ patterns(DF3)      : %.6f (%4d)\n", (double)nb3/DF_SIZE, nb3);
	printf("|\n");
	printf("| < Compact Table 1 Density >\n");
	printf("|   - Compact Table 1(CT1) Size   : %u\n", CT1_TABLE_SIZE);
	printf("|   - Number of non-empty buckets : %d\n", ct1_array_cnt);
	printf("|   - Total num of PIDs in CT1    : %d\n", ct1_pid_cnt);
	printf("|   - Number of PIDs per bucket   : %.3f\n", (double)ct1_pid_cnt/ct1_array_cnt);
	printf("|   - Std Deviation for # of PID  : %.3f\n", ct1_pid_std_dev);
	printf("|\n");
	printf("| < Compact Table 2 Density >\n");
	printf("|   - Compact Table 2(CT2) Size   : %u\n", CT2_TABLE_SIZE);
	printf("|   - Total # of buckets in CT2   : %d\n", ct2_array_cnt);
	printf("|   - Avr # of fragments/bucket   : %.3f (none zero)\n", (double)ct2_array_tot_cnt/ct2_array_cnt);
	printf("|   - Std Dev for # of fragments  : %.3f\n", ct2_array_std_dev);
	printf("|   - Min/Max for # of fragments  : %u / %u\n", ct2_array_min, ct2_array_max);
	printf("|\n");
	printf("|     - Total # of PIDs in CT2    : %d\n", ct2_pid_tot_cnt);
	printf("|     - Avr # of PIDs per frag    : %.3f (none zero)\n", (double)ct2_pid_tot_cnt/ct2_pid_cnt);
	printf("|     - Std Dev for # of PID      : %.3f\n", ct2_pid_std_dev);
	printf("|     - Min/Max for # of PID      : %u / %u\n", ct2_pid_min, ct2_pid_max);
	printf("|\n");
	printf("| < Compact Table 4 Density >\n");
	printf("|   - Compact Table 4(CT4) Size   : %u\n", CT4_TABLE_SIZE);
	printf("|   - Total # of buckets in CT4   : %d\n", ct4_array_cnt);
	printf("|   - Avr # of fragments/bucket   : %.3f (collision rate: %.5f (%d/%d))\n", 
                                                (double)ct4_array_tot_cnt/ct4_array_cnt,
                                                (double)ct4_mt2f_frag_cnt/ct4_array_cnt,
                                                ct4_mt2f_frag_cnt, ct4_array_cnt);
	printf("|   - Std Dev for # of fragments  : %.3f\n", ct4_array_std_dev);
	printf("|   - Min/Max for # of fragments  : %u / %u\n", ct4_array_min, ct4_array_max);
	printf("|\n");
	printf("|     - Total # of PIDs in CT4    : %d\n", ct4_pid_tot_cnt);
	printf("|     - Avr # of PIDs per frag    : %.3f (none zero)\n", (double)ct4_pid_tot_cnt/ct4_pid_cnt);
	printf("|     - Std Dev for # of PID      : %.3f\n", ct4_pid_std_dev);
	printf("|     - Min/Max for # of PID      : %u / %u\n", ct4_pid_min, ct4_pid_max);
	printf("|\n");
	printf("| < Compact Table 8 Density >\n");
	printf("|   - Compact Table 8(CT8) Size   : %u\n", CT8_TABLE_SIZE);
	printf("|   - Total # of buckets in CT8   : %d\n", ct8_array_cnt);
	printf("|   - Avr # of fragments/bucket   : %.3f (collision rate: %.5f (%d/%d))\n", 
                                                (double)ct8_array_tot_cnt/ct8_array_cnt,
                                                (double)ct8_mt2f_frag_cnt/ct8_array_cnt,
                                                ct8_mt2f_frag_cnt, ct8_array_cnt);
	printf("|   - Std Dev for # of fragments  : %.3f\n", ct8_array_std_dev);
	printf("|   - Min/Max for # of fragments  : %u / %u\n", ct8_array_min, ct8_array_max);
	printf("|\n");
	printf("|     - Total # of PIDs in CT8    : %d\n", ct8_pid_tot_cnt);
	printf("|     - Avr # of PIDs per frag    : %.3f (none zero)\n", (double)ct8_pid_tot_cnt/ct8_pid_cnt);
	printf("|     - Std Dev for # of PID      : %.3f\n", ct8_pid_std_dev);
	printf("|     - Min/Max for # of PID      : %u / %u\n", ct8_pid_min, ct8_pid_max);
	printf("|\n");

	if (dfc_total_memory < 1024 * 1024)
	printf("| Total Memory (KB) : %.2f\n", (float)dfc_total_memory/1024);
	else
	printf("| Total Memory (MB) : %.2f\n", (float)dfc_total_memory/(1024*1024));

	if (dfc_pattern_memory < 1024 * 1024)
	printf("|   - Pattern Memory (KB) : %.2f\n", (float)dfc_pattern_memory/1024);
	else
	printf("|   - Pattern Memory (MB) : %.2f\n", (float)dfc_pattern_memory/(1024*1024));

	if (dfc_memory_dfs < 1024 * 1024)
	printf("|   - DF Memory (KB) : %.2f\n", (float)dfc_memory_dfs/1024);
	else
	printf("|   - DF Memory (MB) : %.2f\n", (float)dfc_memory_dfs/(1024*1024));
	
	if (dfc_memory_ct1 < 1024 * 1024)
	printf("|   - CT1 Memory (KB) : %.2f\n", (float)dfc_memory_ct1/1024);
	else
	printf("|   - CT1 Memory (MB) : %.2f\n", (float)dfc_memory_ct1/(1024*1024));
	
	if (dfc_memory_ct2 < 1024 * 1024)
	printf("|   - CT2 Memory (KB) : %.2f\n", (float)dfc_memory_ct2/1024);
	else
	printf("|   - CT2 Memory (MB) : %.2f\n", (float)dfc_memory_ct2/(1024*1024));
	
	if (dfc_memory_ct4 < 1024 * 1024)
	printf("|   - CT4 Memory (KB) : %.2f\n", (float)dfc_memory_ct4/1024);
	else
	printf("|   - CT4 Memory (MB) : %.2f\n", (float)dfc_memory_ct4/(1024*1024));
	
	if (dfc_memory_ct8 < 1024 * 1024)
	printf("|   - CT8 Memory (KB) : %.2f\n", (float)dfc_memory_ct8/1024);
	else
	printf("|   - CT8 Memory (MB) : %.2f\n", (float)dfc_memory_ct8/(1024*1024));
	printf("+--------------------------------------------------------------------------------------\n\n");

	return;

}

static void Add_PID_to_2B_CT(CT_Type_2_2B * CompactTable, uint8_t *temp, PID_TYPE pid, dfcMemoryType type){
    BUC_CNT_TYPE j;
	PID_CNT_TYPE k;
	uint32_t crc = crc64(0, (uint64_t)temp, 2); //_mm_crc32_u16(0, *(uint16_t*)temp);

	crc &= CT2_TABLE_SIZE_MASK;

	if( CompactTable[crc].cnt != 0){
		for(j = 0; j < CompactTable[crc].cnt; j++){
			if( CompactTable[crc].array[j].pat == *(uint16_t*)temp )
				break;
		}

		if( j == CompactTable[crc].cnt){ // If not found,
			CompactTable[crc].cnt++;
			CompactTable[crc].array =
				(CT_Type_2_2B_Array *)DFC_REALLOC((void*)CompactTable[crc].array,
							CompactTable[crc].cnt, DFC_CT_Type_2_2B_Array, type);
			CompactTable[crc].array[CompactTable[crc].cnt-1].pat = *(uint16_t*)temp;

			CompactTable[crc].array[CompactTable[crc].cnt-1].cnt = 1;
			CompactTable[crc].array[CompactTable[crc].cnt-1].pid =
									(PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE), type);
			CompactTable[crc].array[CompactTable[crc].cnt-1].pid[0] = pid;
		}else{ // If found,
			for(k = 0; k < CompactTable[crc].array[j].cnt; k++){
				if(CompactTable[crc].array[j].pid[k] == pid)
					break;
			}
			if( k == CompactTable[crc].array[j].cnt ){
				CompactTable[crc].array[j].cnt++;
				CompactTable[crc].array[j].pid =
					 (PID_TYPE *)DFC_REALLOC((void*)CompactTable[crc].array[j].pid,
						CompactTable[crc].array[j].cnt, DFC_PID_TYPE, type);
				CompactTable[crc].array[j].pid[CompactTable[crc].array[j].cnt-1] =
					 pid;
			}
		}
	}else{ // If there is no elements in the CT4,
		CompactTable[crc].cnt = 1;
		CompactTable[crc].array = (CT_Type_2_2B_Array *)DFC_MALLOC(sizeof(CT_Type_2_2B_Array), type);
		memset(CompactTable[crc].array, 0, sizeof(CT_Type_2_2B_Array));

		CompactTable[crc].array[0].pat = *(uint16_t*)temp;
		CompactTable[crc].array[0].cnt = 1;
		CompactTable[crc].array[0].pid = (PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE), type);
		CompactTable[crc].array[0].pid[0] = pid;
	}	

}

int DFC_Compile(DFC_STRUCTURE* dfc)
{
	uint32_t i=0;
    uint32_t alpha_cnt;

    int j, k, l;
    BUC_CNT_TYPE m, n;
	DFC_PATTERN *plist;

    uint8_t temp[8], flag[8];
    uint16_t fragment_16;
    uint32_t fragment_32;
    uint64_t fragment_64;
    uint32_t byteIndex, bitMask;

	dfc_memory_ct1 = sizeof(CT_Type_1)    * CT1_TABLE_SIZE;
	dfc_memory_ct2 = 0;
	dfc_memory_ct3 = 0;
	dfc_memory_ct4 = 0;
	dfc_memory_ct8 = 0;
	dfc_total_memory = sizeof(DFC_STRUCTURE) + dfc_pattern_memory;

/* ####################################################################################### */
/* ###############                  MatchList initialization              ################ */
/* ####################################################################################### */
    int begin_node_flag = 1;
    for (i = 0; i < INIT_HASH_SIZE; i++) {
        DFC_PATTERN *node = dfc->init_hash[i], *prev_node;
        int first_node_flag = 1;
        while(node != NULL) {
            if(begin_node_flag){
                begin_node_flag = 0;
                dfc->dfcPatterns = node;
            }else{
                if(first_node_flag){
                    first_node_flag = 0;
                    prev_node->next = node;
                }
            }
            prev_node = node;
            node = node->next;
        }
    }

    free(dfc->init_hash);
    dfc->init_hash = NULL;

	dfc->dfcMatchList = (DFC_PATTERN **)DFC_MALLOC(sizeof(DFC_PATTERN*) * dfc->numPatterns, 
									    						DFC_MEMORY_TYPE__PATTERN);
    MEMASSERT_DFC(dfc->dfcMatchList, "DFC_Compile");

	for(plist = dfc->dfcPatterns; plist != NULL; plist = plist->next){
        if(dfc->dfcMatchList[plist->iid] != NULL){
            printe("Internal ID ERROR : %u\n", plist->iid);
        }
		dfc->dfcMatchList[plist->iid] = plist;
	}


/* ####################################################################################### */

/* ####################################################################################### */
/* ###############              0. Direct Filters initialization          ################ */
/* ####################################################################################### */

    /* Initializing Bloom Filter */
    for (i = 0; i < DF_SIZE_REAL; i++) {
        dfc->DirectFilter1[i] = 0;
        dfc->ADD_DF_4_plus[i] = 0;
        dfc->ADD_DF_8_1[i] = 0;
        dfc->ADD_DF_4_1[i] = 0;
        dfc->cDF2[i] = 0;
        dfc->ADD_DF_8_2[i] = 0;
        dfc->cDF1[i] = 0;
    }

    for(i = 0; i < 256; i++){
        dfc->cDF0[i] = 0;
    }

/* ####################################################################################### */

/* ####################################################################################### */
/* ###############            1. Decide extracting position               ################ */
/* ####################################################################################### */
	pattern_interval = 32;
    min_pattern_interval = 0;
#if 0
    int min_bit_count = 65536;
    int max_bucket_cnt = 0;

    uint8_t tmpADD_DF_8_1_2[DF_SIZE_REAL];
    uint8_t tmpADD_DF_8_1_3[DF_SIZE_REAL];
    uint8_t tmpADD_DF_8_1_4[DF_SIZE_REAL];
    uint8_t tmpDirectFilter1[DF_SIZE_REAL];
	CT_Type_2_8B tmpCompactTable8[CT8_TABLE_SIZE];

    for(i = 0; i <= pattern_interval; i++)
    {
        int l;
        int pcount = 0;
		memset(tmpDirectFilter1, 0, sizeof(uint8_t) * DF_SIZE_REAL);
		memset(tmpADD_DF_8_1_2, 0, sizeof(uint8_t) * DF_SIZE_REAL);
		memset(tmpADD_DF_8_1_3, 0, sizeof(uint8_t) * DF_SIZE_REAL);
		memset(tmpADD_DF_8_1_4, 0, sizeof(uint8_t) * DF_SIZE_REAL);
		memset(tmpCompactTable8, 0, sizeof(CT_Type_2_8B) * CT8_TABLE_SIZE);

        for(plist = dfc->dfcPatterns; plist != NULL; plist = plist->next){
            if (plist->n > 1) {
                alpha_cnt = 0;

                do {
                    for (j=1, k=0; j>=0; --j, k++) {
                        flag[k] =  (alpha_cnt >> j) & 1;
                    }

                    if (plist->n == 2) {
                        for (j=plist->n - 2, k=0; j< plist->n; j++, k++){
                            Build_pattern(plist, flag, temp, i,j,k);
                        }
                    }else if (plist->n == 3) {
                        for (j=0 , k=0; j < 2; j++, k++){
                            Build_pattern(plist, flag, temp, i,j,k);
                        }
                    }else if (plist->n < 8) {
                        for (j=plist->n - 4, k=0; j< plist->n-2; j++, k++){
                            Build_pattern(plist, flag, temp, i,j,k);
                        }
                    } else { // len >= 8
						uint8_t temp2[8];
                        for (j = i*(plist->n-8)/pattern_interval, k=0; 
							 j < i*(plist->n-8)/pattern_interval+8; j++,k++){
                            Build_pattern(plist, flag, temp, i,j,k);
							temp2[k] = plist->patrn[j];
                        }
						uint64_t crc = _mm_crc32_u64(0, *(uint64_t*)temp2);
						crc &= CT8_TABLE_SIZE_MASK;

#if 1
						if( tmpCompactTable8[crc].cnt != 0){ 
							for(j = 0; j < tmpCompactTable8[crc].cnt; j++){
								if( tmpCompactTable8[crc].array[j].pat == *(uint64_t*)temp2 )
									break;	
							}

							if( j == tmpCompactTable8[crc].cnt){ // If not found,
								tmpCompactTable8[crc].cnt++;
								tmpCompactTable8[crc].array = 
									(CT_Type_2_8B_Array *)DFC_REALLOC((void*)tmpCompactTable8[crc].array,
									tmpCompactTable8[crc].cnt, DFC_CT_Type_2_8B_Array, DFC_MEMORY_TYPE__NONE);
								tmpCompactTable8[crc].
										array[tmpCompactTable8[crc].cnt-1].pat=*(uint64_t*)temp2;

								tmpCompactTable8[crc].array[tmpCompactTable8[crc].cnt-1].cnt = 1;
							}
						}else{ // If there is no elements in the CT8,
							tmpCompactTable8[crc].cnt = 1;
							tmpCompactTable8[crc].array =
									 (CT_Type_2_8B_Array *)DFC_MALLOC(sizeof(CT_Type_2_8B_Array),
																				 DFC_MEMORY_TYPE__NONE);
							memset(tmpCompactTable8[crc].array, 0, sizeof(CT_Type_2_8B_Array));

							tmpCompactTable8[crc].array[0].pat = *(uint64_t*)temp2;
							tmpCompactTable8[crc].array[0].cnt = 1;
						} 
#endif
                    }

                    byteIndex = (uint32_t)BINDEX((*(uint16_t*)temp) & DF_MASK);
                    bitMask = BMASK((*(uint16_t*)temp) & DF_MASK);

                    tmpDirectFilter1[byteIndex] |= bitMask;
                    if(!i)
                    {
                        if (plist->n == 2) {
                            dfc->DirectFilter5[byteIndex] |= bitMask;
                        } else if (plist->n == 3) {
                            dfc->DirectFilter6[byteIndex] |= bitMask;
                        }
                    }
                    alpha_cnt++;
                } while (alpha_cnt < 4);
            }
            if (plist->n >= 8) {
                alpha_cnt = 0;
                pcount++;
                do {
                    for (j=7, k=0; j>=0; --j, k++) {
                        flag[k] =  (alpha_cnt >> j) & 1;
                    }

                    for (j=i*(plist->n - 8)/pattern_interval, k=0; j< i*(plist->n-8)/pattern_interval+8; j++, k++){
                        Build_pattern(plist, flag, temp, i,j,k);
                    }

                    byteIndex = BINDEX((*(((uint16_t*)temp)+1)) & DF_MASK);
                    bitMask = BMASK((*(((uint16_t*)temp)+1)) & DF_MASK);
                    tmpADD_DF_8_1_2[byteIndex] |= bitMask;

                    byteIndex = BINDEX((*(((uint16_t*)temp)+2)) & DF_MASK);
                    bitMask = BMASK((*(((uint16_t*)temp)+2)) & DF_MASK);
                    tmpADD_DF_8_1_3[byteIndex] |= bitMask;

                    byteIndex = BINDEX((*(((uint16_t*)temp)+3)) & DF_MASK);
                    bitMask = BMASK((*(((uint16_t*)temp)+3)) & DF_MASK);
                    tmpADD_DF_8_1_4[byteIndex] |= bitMask;

                    alpha_cnt++;
                } while (alpha_cnt < 256);
            }
        }
        int bit_count = 0;
        int bit_count_8_2 = 0;
        int bit_count_8_3 = 0;
        int bit_count_8_4 = 0;
        int bit_count_ct8 = 0;
        for(l = 0; l < DF_SIZE_REAL; l++)
        {
            bit_count += __builtin_popcount(tmpDirectFilter1[l]);
            bit_count_8_2 += __builtin_popcount(tmpADD_DF_8_1_2[l]);
            bit_count_8_3 += __builtin_popcount(tmpADD_DF_8_1_3[l]);
            bit_count_8_4 += __builtin_popcount(tmpADD_DF_8_1_4[l]);
        }
        for(l = 0; l < CT8_TABLE_SIZE; l++)
        {
            //bit_count_ct8 += __builtin_popcount(tmpCompactTable8[l]);
            bit_count_ct8 += tmpCompactTable8[l].cnt;
			DFC_FREE(tmpCompactTable8[l].array, sizeof(CT_Type_2_8B_Array) * tmpCompactTable8[l].cnt);
		}
        //if(bit_count < min_bit_count && (double)bit_count_8/pcount < 1.45)
        //if(bit_count < min_bit_count)
        //if(bit_count_ct8 < max_bucket_cnt)
        if(bit_count_8_3 + bit_count_8_4 < min_bit_count)
        {
            min_bit_count = bit_count_8_3 + bit_count_8_4;
            min_pattern_interval = i;
        }
		//printf("position: %d, DF1: %d, DF3_2 (8B): %d, DF3_3 (8B): %d, DF3_4 (8B): %d # of buckets in CT8: %d\n",
		//			 i, bit_count, bit_count_8_2, bit_count_8_3, bit_count_8_4 , bit_count_ct8);
    }
#endif
	min_pattern_interval = 32;
	//printf("Extracting Position : %d/%d\n", min_pattern_interval, pattern_interval);



/* ####################################################################################### */

/* ####################################################################################### */
/* ###############               Direct Filters setup                     ################ */
/* ####################################################################################### */
    memset(dfc->CompactTable1, 0, sizeof(CT_Type_1)*CT1_TABLE_SIZE);

	for(plist = dfc->dfcPatterns; plist != NULL; plist = plist->next){
		/* 0. Initialization for DF8 (for 1B patterns)*/
        if (plist->n == 1) {
            temp[0] = plist->casepatrn[0];
            for(j = 0; j < 256; j++){
                temp[1] = j;

                fragment_16 = (temp[1] << 8) | temp[0];
                byteIndex = (uint32_t)BINDEX(fragment_16 & DF_MASK);
                bitMask = BMASK(fragment_16 & DF_MASK);

                dfc->DirectFilter1[byteIndex] |= bitMask;
            }

            dfc->cDF0[temp[0]] = 1;
			if(dfc->CompactTable1[temp[0]].cnt == 0){
				dfc->CompactTable1[temp[0]].cnt++;
				dfc->CompactTable1[temp[0]].pid[0] = plist->iid;
			}else{
				for(k = 0; k < dfc->CompactTable1[temp[0]].cnt; k++){
					if(dfc->CompactTable1[temp[0]].pid[k] == plist->iid)
						break;
				}
				if(k == dfc->CompactTable1[temp[0]].cnt){
					dfc->CompactTable1[temp[0]].pid[dfc->CompactTable1[temp[0]].cnt++] = plist->iid;
					if(dfc->CompactTable1[temp[0]].cnt >= CT_TYPE1_PID_CNT_MAX)
						printf("Too many PIDs in CT1. You should expand the size.\n");
				}
			}

			if(plist->nocase){
				if(plist->casepatrn[0] >= 97/*a*/ && plist->casepatrn[0] <= 122/*z*/){
					/* when the pattern is lower case */
					temp[0] = toupper(plist->casepatrn[0]);
				}else{	
					/* when the pattern is upper case */
					temp[0] = tolower(plist->casepatrn[0]);
				}

				for(j = 0; j < 256; j++){
					temp[1] = j;

                    fragment_16 = (temp[1] << 8) | temp[0];
					byteIndex = (uint32_t)BINDEX(fragment_16 & DF_MASK);
					bitMask = BMASK(fragment_16 & DF_MASK);

					dfc->DirectFilter1[byteIndex] |= bitMask;
				}

				dfc->cDF0[temp[0]] = 1;
				if(dfc->CompactTable1[temp[0]].cnt == 0){
					dfc->CompactTable1[temp[0]].cnt++;
					dfc->CompactTable1[temp[0]].pid[0] = plist->iid;
				}else{
					for(k = 0; k < dfc->CompactTable1[temp[0]].cnt; k++){
						if(dfc->CompactTable1[temp[0]].pid[k] == plist->iid)
							break;
					}
					if(k == dfc->CompactTable1[temp[0]].cnt){
						dfc->CompactTable1[temp[0]].pid[dfc->CompactTable1[temp[0]].cnt++] = plist->iid;
						if(dfc->CompactTable1[temp[0]].cnt >= CT_TYPE1_PID_CNT_MAX)
							printf("Too many PIDs in CT1. You should expand the size.\n");
					}
				}
			}
        }

	    /* 1. Initialization for DF1 */
        if (plist->n > 1) {
            alpha_cnt = 0;

            do {
                for (j=1, k=0; j>=0; --j, k++) {
                    flag[k] =  (alpha_cnt >> j) & 1;
                }

				if (plist->n == 2) {
                    for (j=plist->n - 2, k=0; j< plist->n; j++, k++){
                        Build_pattern(plist, flag, temp, i,j,k);
                    }
				}else if (plist->n == 3) {
                    //for (j=0 , k=0; j < 2; j++, k++){
                    for (j=plist->n - 2, k=0; j< plist->n; j++, k++){
                        Build_pattern(plist, flag, temp, i,j,k);
                    }
				}else if (plist->n < 8) {
                    for (j=plist->n - 4, k=0; j< plist->n-2; j++, k++){
                        Build_pattern(plist, flag, temp, i,j,k);
                    }
                } else { // len >= 8
                    for (j = min_pattern_interval*(plist->n-8)/pattern_interval, k = 0; 
						 j < min_pattern_interval*(plist->n - 8)/pattern_interval + 2; j++, k++){
                        Build_pattern(plist, flag, temp, i,j,k);
                    }
                }

                fragment_16 = (temp[1] << 8) | temp[0];
                byteIndex = (uint32_t)BINDEX(fragment_16 & DF_MASK);
                bitMask = BMASK(fragment_16 & DF_MASK);

                dfc->DirectFilter1[byteIndex] |= bitMask;

				if(plist->n == 2 || plist->n == 3)
					dfc->cDF1[byteIndex] |= bitMask;

                alpha_cnt++;
            } while (alpha_cnt < 4);
        }

		/* Initializing 4B DF, 8B DF */
        if (plist->n >= 4) {
            alpha_cnt = 0;

            do {
                for (j=3, k=0; j>=0; --j, k++) {
                    flag[k] =  (alpha_cnt >> j) & 1;
                }

                if (plist->n < 8) {
                    for (j=plist->n - 4, k=0; j< plist->n; j++, k++) {
                        Build_pattern(plist, flag, temp, i,j,k);
                    }
                } else {
                    for (j = min_pattern_interval*(plist->n-8)/pattern_interval, k=0;
						 j < min_pattern_interval*(plist->n-8)/pattern_interval + 4; j++, k++) {
                        Build_pattern(plist, flag, temp, i,j,k);
                    }
                }

                byteIndex = BINDEX((*(((uint16_t*)temp)+1)) & DF_MASK);
                bitMask = BMASK((*(((uint16_t*)temp)+1)) & DF_MASK);

                dfc->ADD_DF_4_plus[byteIndex] |= bitMask;
                if (plist->n >= 4 && plist->n < 8) {
                    dfc->ADD_DF_4_1[byteIndex] |= bitMask;

                    fragment_16 = (temp[1] << 8) | temp[0];
                    byteIndex = BINDEX(fragment_16 & DF_MASK);
                    bitMask = BMASK(fragment_16 & DF_MASK);

                    dfc->cDF2[byteIndex] |= bitMask;
				}
                alpha_cnt++;
            } while (alpha_cnt < 16);
        }

        if (plist->n >= 8) {
            alpha_cnt = 0;

            do {
                for (j=7, k=0; j>=0; --j, k++) {
                    flag[k] =  (alpha_cnt >> j) & 1;
                }

                for (j = min_pattern_interval*(plist->n - 8)/pattern_interval, k=0; 
					 j < min_pattern_interval*(plist->n - 8)/pattern_interval + 8; j++, k++){
                    Build_pattern(plist, flag, temp, i,j,k);
                }

                byteIndex = BINDEX((*(((uint16_t*)temp)+3)) & DF_MASK);
                bitMask = BMASK((*(((uint16_t*)temp)+3)) & DF_MASK);

                dfc->ADD_DF_8_1[byteIndex] |= bitMask;

                byteIndex = BINDEX((*(((uint16_t*)temp)+2)) & DF_MASK);
                bitMask = BMASK((*(((uint16_t*)temp)+2)) & DF_MASK);

                dfc->ADD_DF_8_2[byteIndex] |= bitMask;

                alpha_cnt++;
            } while (alpha_cnt < 256);
        }

    }

    //printf("DF Initialization is done.\n");
/* ####################################################################################### */

/* ####################################################################################### */
/* ###############                Compact Tables initialization           ################ */
/* ####################################################################################### */

	dfc_memory_ct2 += sizeof(CT_Type_2)*CT2_TABLE_SIZE;
    memset(dfc->CompactTable2, 0, sizeof(CT_Type_2)*CT2_TABLE_SIZE);

	dfc_memory_ct4 += sizeof(CT_Type_2)*CT4_TABLE_SIZE;
    memset(dfc->CompactTable4, 0, sizeof(CT_Type_2)*CT4_TABLE_SIZE);

	dfc_memory_ct8 += sizeof(CT_Type_2_8B)*CT8_TABLE_SIZE;
    memset(dfc->CompactTable8, 0, sizeof(CT_Type_2_8B)*CT8_TABLE_SIZE);


/* ####################################################################################### */

/* ####################################################################################### */
/* ###############                   Compact Tables setup                 ################ */
/* ####################################################################################### */

	for(plist = dfc->dfcPatterns; plist != NULL; plist = plist->next){
        if (plist->n == 2 || plist->n == 3) {
			alpha_cnt = 0;

            do {
                for (j=1, k=0; j>=0; --j, k++) {
                    flag[k] =  (alpha_cnt >> j) & 1;
                }

                for (j=plist->n - 2, k=0; j< plist->n; j++, k++){
                    Build_pattern(plist, flag, temp, i,j,k);
                }

				 // 2.
                fragment_16 = (temp[1] << 8) | temp[0];
                uint32_t crc = crc64(0, (uint64_t)fragment_16, 2); //_mm_crc32_u16(0, fragment_16);

                // 3.
                crc &= CT2_TABLE_SIZE_MASK;

                // 4.
                if( dfc->CompactTable2[crc].cnt != 0){
                    for(n = 0; n < dfc->CompactTable2[crc].cnt; n++){
                        if( dfc->CompactTable2[crc].array[n].pat == fragment_16 )
                            break;
                    }

                    if( n == dfc->CompactTable2[crc].cnt){ // If not found,
                        dfc->CompactTable2[crc].cnt++;
                        dfc->CompactTable2[crc].array =
                            (CT_Type_2_Array *)DFC_REALLOC((void*)dfc->CompactTable2[crc].array,
                                        dfc->CompactTable2[crc].cnt, DFC_CT_Type_2_Array, DFC_MEMORY_TYPE__CT2);
                        dfc->CompactTable2[crc].array[dfc->CompactTable2[crc].cnt-1].pat = fragment_16;

                        dfc->CompactTable2[crc].array[dfc->CompactTable2[crc].cnt-1].cnt = 1;
                        dfc->CompactTable2[crc].array[dfc->CompactTable2[crc].cnt-1].pid =
                                                (PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE), DFC_MEMORY_TYPE__CT2);
                        dfc->CompactTable2[crc].array[dfc->CompactTable2[crc].cnt-1].pid[0] = plist->iid;
						dfc->CompactTable2[crc].array[dfc->CompactTable2[crc].cnt-1].DirectFilter = NULL;
						dfc->CompactTable2[crc].array[dfc->CompactTable2[crc].cnt-1].CompactTable = NULL;
                    }else{ // If found,
                        for(m = 0; m < dfc->CompactTable2[crc].array[n].cnt; m++){
                            if(dfc->CompactTable2[crc].array[n].pid[m] == plist->iid)
                                break;
                        }
                        if( m == dfc->CompactTable2[crc].array[n].cnt ){
                            dfc->CompactTable2[crc].array[n].cnt++;
                            dfc->CompactTable2[crc].array[n].pid =
                                 (PID_TYPE *)DFC_REALLOC((void*)dfc->CompactTable2[crc].array[n].pid,
                                    dfc->CompactTable2[crc].array[n].cnt, DFC_PID_TYPE, DFC_MEMORY_TYPE__CT2);
                            dfc->CompactTable2[crc].array[n].pid[dfc->CompactTable2[crc].array[n].cnt-1] =
                                 plist->iid;
                        }
                    }
                }else{ // If there is no elements in the CT4,
                    dfc->CompactTable2[crc].cnt = 1;
                    dfc->CompactTable2[crc].array = (CT_Type_2_Array *)DFC_MALLOC(sizeof(CT_Type_2_Array),
                                                                                             DFC_MEMORY_TYPE__CT2);
                    memset(dfc->CompactTable2[crc].array, 0, sizeof(CT_Type_2_Array));

                    dfc->CompactTable2[crc].array[0].pat = fragment_16;
                    dfc->CompactTable2[crc].array[0].cnt = 1;
                    dfc->CompactTable2[crc].array[0].pid = (PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE),
                                                                             DFC_MEMORY_TYPE__CT2);
                    dfc->CompactTable2[crc].array[0].pid[0] = plist->iid;
					dfc->CompactTable2[crc].array[0].DirectFilter = NULL;
					dfc->CompactTable2[crc].array[0].CompactTable = NULL;
                }

                alpha_cnt++;
            } while (alpha_cnt < 4);
        }

		/* CT4 initialization */
        if (plist->n >= 4 && plist->n < 8) {
            alpha_cnt = 0;
            do {
				// 1. 

				for (j=3, k=0; j>=0; --j, k++) {
					flag[k] =  (alpha_cnt >> j) & 1;
				}

				for (j=plist->n - 4, k=0; j< plist->n; j++, k++){
					Build_pattern(plist, flag, temp, i,j,k);
				}

				// 2.		
                fragment_32 = (temp[3] << 24) | (temp[2] << 16) | (temp[1] << 8) | temp[0];
				uint32_t crc = crc64(0, (uint64_t)fragment_32, 4); //_mm_crc32_u32(0, fragment_32);

				// 3.
				crc &= CT4_TABLE_SIZE_MASK;			

				// 4.
				if( dfc->CompactTable4[crc].cnt != 0){ 
					for(n = 0; n < dfc->CompactTable4[crc].cnt; n++){
						if( dfc->CompactTable4[crc].array[n].pat == fragment_32 )
							break;	
					}

					if( n == dfc->CompactTable4[crc].cnt){ // If not found,
						dfc->CompactTable4[crc].cnt++;
						dfc->CompactTable4[crc].array = 
							(CT_Type_2_Array *)DFC_REALLOC((void*)dfc->CompactTable4[crc].array,
										dfc->CompactTable4[crc].cnt, DFC_CT_Type_2_Array, DFC_MEMORY_TYPE__CT4);
						dfc->CompactTable4[crc].array[dfc->CompactTable4[crc].cnt-1].pat = fragment_32;

						dfc->CompactTable4[crc].array[dfc->CompactTable4[crc].cnt-1].cnt = 1;
						dfc->CompactTable4[crc].array[dfc->CompactTable4[crc].cnt-1].pid = 
												(PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE), DFC_MEMORY_TYPE__CT4);
						dfc->CompactTable4[crc].array[dfc->CompactTable4[crc].cnt-1].pid[0] = plist->iid;
						dfc->CompactTable4[crc].array[dfc->CompactTable4[crc].cnt-1].DirectFilter = NULL;
						dfc->CompactTable4[crc].array[dfc->CompactTable4[crc].cnt-1].CompactTable = NULL;
					}else{ // If found,
						for(m = 0; m < dfc->CompactTable4[crc].array[n].cnt; m++){
							if(dfc->CompactTable4[crc].array[n].pid[m] == plist->iid)
								break;
						}
						if( m == dfc->CompactTable4[crc].array[n].cnt ){
							dfc->CompactTable4[crc].array[n].cnt++;
							dfc->CompactTable4[crc].array[n].pid =
								 (PID_TYPE *)DFC_REALLOC((void*)dfc->CompactTable4[crc].array[n].pid,
									dfc->CompactTable4[crc].array[n].cnt, DFC_PID_TYPE, DFC_MEMORY_TYPE__CT4);
							dfc->CompactTable4[crc].array[n].pid[dfc->CompactTable4[crc].array[n].cnt-1] =
								 plist->iid;
						}
					}
				}else{ // If there is no elements in the CT4,
					dfc->CompactTable4[crc].cnt = 1;
					dfc->CompactTable4[crc].array = (CT_Type_2_Array *)DFC_MALLOC(sizeof(CT_Type_2_Array),
																							 DFC_MEMORY_TYPE__CT4);
					memset(dfc->CompactTable4[crc].array, 0, sizeof(CT_Type_2_Array));

					dfc->CompactTable4[crc].array[0].pat = fragment_32;
					dfc->CompactTable4[crc].array[0].cnt = 1;
					dfc->CompactTable4[crc].array[0].pid = (PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE),
																			 DFC_MEMORY_TYPE__CT4);
					dfc->CompactTable4[crc].array[0].pid[0] = plist->iid;
					dfc->CompactTable4[crc].array[0].DirectFilter = NULL;
					dfc->CompactTable4[crc].array[0].CompactTable = NULL;
				} 
                alpha_cnt++;
            } while (alpha_cnt < 16);
        }

		/* CT8 initialization */
        if (plist->n >= 8) {
            alpha_cnt = 0;
            do {
                for (j=7, k=0; j>=0; --j, k++) {
                    flag[k] =  (alpha_cnt >> j) & 1;
                }


#ifdef CT8_SWITCH
				for (j = min_pattern_interval*(plist->n - 8)/pattern_interval, k=0;
					 j < min_pattern_interval*(plist->n - 8)/pattern_interval + 8; j++, k++) {
					temp[k] = plist->patrn[j];
				}
#else
                for (j=plist->n - 8, k=0; j< plist->n; j++, k++){
                    Build_pattern(plist, flag, temp, i,j,k);
                }
#endif

                // 1. Calulating Indice
                fragment_32 = (temp[7] << 24) | (temp[6] << 16) | (temp[5] << 8) | temp[4];
                fragment_64 = ((uint64_t)fragment_32 << 32) |
                              (temp[3] << 24) | (temp[2] << 16) | (temp[1] << 8) | temp[0];
				uint64_t crc = crc64(0, (uint64_t)fragment_64, 8); //_mm_crc32_u64(0, fragment_64);
				crc &= CT8_TABLE_SIZE_MASK;

				if( dfc->CompactTable8[crc].cnt != 0){ 
					for(n = 0; n < dfc->CompactTable8[crc].cnt; n++){
						if( dfc->CompactTable8[crc].array[n].pat == fragment_64)
							break;	
					}

					if( n == dfc->CompactTable8[crc].cnt){ // If not found,
						dfc->CompactTable8[crc].cnt++;
						dfc->CompactTable8[crc].array = 
							(CT_Type_2_8B_Array *)DFC_REALLOC((void*)dfc->CompactTable8[crc].array,
										dfc->CompactTable8[crc].cnt, DFC_CT_Type_2_8B_Array, DFC_MEMORY_TYPE__CT8);
						dfc->CompactTable8[crc].array[dfc->CompactTable8[crc].cnt-1].pat = fragment_64;

						dfc->CompactTable8[crc].array[dfc->CompactTable8[crc].cnt-1].cnt = 1;
						dfc->CompactTable8[crc].array[dfc->CompactTable8[crc].cnt-1].pid = 
												(PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE), DFC_MEMORY_TYPE__CT8);
						dfc->CompactTable8[crc].array[dfc->CompactTable8[crc].cnt-1].pid[0] = plist->iid;
						dfc->CompactTable8[crc].array[dfc->CompactTable8[crc].cnt-1].DirectFilter = NULL;
						dfc->CompactTable8[crc].array[dfc->CompactTable8[crc].cnt-1].CompactTable = NULL;
					}else{ // If found,
						for(m = 0; m < dfc->CompactTable8[crc].array[n].cnt; m++){
							if(dfc->CompactTable8[crc].array[n].pid[m] == plist->iid)
								break;
						}
						if( m == dfc->CompactTable8[crc].array[n].cnt ){
							dfc->CompactTable8[crc].array[n].cnt++;
							dfc->CompactTable8[crc].array[n].pid =
								 (PID_TYPE *)DFC_REALLOC((void*)dfc->CompactTable8[crc].array[n].pid,
									dfc->CompactTable8[crc].array[n].cnt, DFC_PID_TYPE, DFC_MEMORY_TYPE__CT8);
							dfc->CompactTable8[crc].array[n].pid[dfc->CompactTable8[crc].array[n].cnt-1] =
																									 plist->iid;
						}
					}
				}else{ // If there is no elements in the CT8,
					dfc->CompactTable8[crc].cnt = 1;
					dfc->CompactTable8[crc].array = (CT_Type_2_8B_Array *)DFC_MALLOC(sizeof(CT_Type_2_8B_Array),
																							 DFC_MEMORY_TYPE__CT8);
					memset(dfc->CompactTable8[crc].array, 0, sizeof(CT_Type_2_8B_Array));

					dfc->CompactTable8[crc].array[0].pat = fragment_64;
					dfc->CompactTable8[crc].array[0].cnt = 1;
					dfc->CompactTable8[crc].array[0].pid = (PID_TYPE *)DFC_MALLOC(sizeof(PID_TYPE),
																			 DFC_MEMORY_TYPE__CT8);
					dfc->CompactTable8[crc].array[0].pid[0] = plist->iid;
					dfc->CompactTable8[crc].array[0].DirectFilter = NULL;
					dfc->CompactTable8[crc].array[0].CompactTable = NULL;
				} 
                alpha_cnt++;
            } while (alpha_cnt < 256);
        }
    }
    //printf("CT Initialization is done.\n");
/* ####################################################################################### */

/* ####################################################################################### */
/* ###############                   Recursive filtering                  ################ */
/* ####################################################################################### */
#ifdef ENABLE_RECURSIVE
	// Only for CT2 firstly
#if 1
	for(i=0; i < CT2_TABLE_SIZE; i++){
		for(n=0; n < dfc->CompactTable2[i].cnt; n++){
			/* If the number of PID is bigger than 3, do recursive filtering */
			if(dfc->CompactTable2[i].array[n].cnt >= RECURSIVE_BOUNDARY){
				/* Initialization */
				dfc->CompactTable2[i].array[n].DirectFilter=(uint8_t*)DFC_MALLOC(sizeof(uint8_t)*DF_SIZE_REAL,
																						DFC_MEMORY_TYPE__CT2);
				dfc->CompactTable2[i].array[n].CompactTable=(CT_Type_2_2B*)DFC_MALLOC
                                            (sizeof(CT_Type_2_2B)*RECURSIVE_CT_SIZE, DFC_MEMORY_TYPE__CT2);

				if(!dfc->CompactTable2[i].array[n].DirectFilter || !dfc->CompactTable2[i].array[n].CompactTable){
					printe("Failed to allocate memory for recursive things.\n");
					return -1;
				}

				PID_TYPE *tempPID = (PID_TYPE*)DFC_MALLOC(sizeof(PID_TYPE)*dfc->CompactTable2[i].array[n].cnt,
																					    DFC_MEMORY_TYPE__CT2);	
				memcpy(tempPID, dfc->CompactTable2[i].array[n].pid, 
														sizeof(PID_TYPE) * dfc->CompactTable2[i].array[n].cnt);
				//free(dfc->CompactTable2[i].array[n].pid);
				DFC_FREE(dfc->CompactTable2[i].array[n].pid, dfc->CompactTable2[i].array[n].cnt*sizeof(PID_TYPE),
															DFC_MEMORY_TYPE__CT2);
				dfc->CompactTable2[i].array[n].pid = NULL;

				int temp_cnt = 0; // cnt for 2 byte patterns.

				for(m = 0; m < dfc->CompactTable2[i].array[n].cnt; m++){
					//DFC_PATTERN *mlist = dfc->dfcMatchList[]; 
					int pat_len = dfc->dfcMatchList[tempPID[m]]->n - 2;

					if(pat_len == 0){ /* When pat length is 2 */
						temp_cnt ++;
						dfc->CompactTable2[i].array[n].pid =
							 (PID_TYPE *)realloc(dfc->CompactTable2[i].array[n].pid, sizeof(PID_TYPE)*temp_cnt);
						dfc->CompactTable2[i].array[n].pid[temp_cnt-1] = tempPID[m];
					}else if (pat_len == 1){ /* When pat length is 3 */
						if(dfc->dfcMatchList[tempPID[m]]->nocase){
							temp[1] = tolower(dfc->dfcMatchList[tempPID[m]]->patrn[0]);
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
								byteIndex = BINDEX(fragment_16);
								bitMask = BMASK(fragment_16);

								dfc->CompactTable2[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable2[i].array[n].CompactTable, temp,	tempPID[m],
																					DFC_MEMORY_TYPE__CT2);
							}

							temp[1] = toupper(dfc->dfcMatchList[tempPID[m]]->patrn[0]);
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
								byteIndex = BINDEX(fragment_16);
								bitMask = BMASK(fragment_16);

								dfc->CompactTable2[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable2[i].array[n].CompactTable, temp,	tempPID[m],
																					DFC_MEMORY_TYPE__CT2);
							}
						}else{
							temp[1] = dfc->dfcMatchList[tempPID[m]]->casepatrn[0];
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable2[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable2[i].array[n].CompactTable, temp,	tempPID[m],
																					DFC_MEMORY_TYPE__CT2);
							}
						}
					}
				}

				dfc->CompactTable2[i].array[n].cnt = temp_cnt;
				DFC_FREE(tempPID, sizeof(PID_TYPE)*dfc->CompactTable2[i].array[n].cnt, DFC_MEMORY_TYPE__CT2);
			}
		}
	}
#endif
	// Only for CT4 firstly
	for(i=0; i < CT4_TABLE_SIZE; i++){
		for(n=0; n < dfc->CompactTable4[i].cnt; n++){
			/* If the number of PID is bigger than 3, do recursive filtering */
			if(dfc->CompactTable4[i].array[n].cnt >= RECURSIVE_BOUNDARY){
				/* Initialization */
				dfc->CompactTable4[i].array[n].DirectFilter=(uint8_t*)DFC_MALLOC(sizeof(uint8_t)*DF_SIZE_REAL,
																							DFC_MEMORY_TYPE__CT4);
				dfc->CompactTable4[i].array[n].CompactTable=(CT_Type_2_2B*)DFC_MALLOC
                                                (sizeof(CT_Type_2_2B)*RECURSIVE_CT_SIZE, DFC_MEMORY_TYPE__CT4);

				if(!dfc->CompactTable4[i].array[n].DirectFilter || !dfc->CompactTable4[i].array[n].CompactTable){
					printe("Failed to allocate memory for recursive things.\n");
					return -1;
				}

				PID_TYPE *tempPID = (PID_TYPE*)DFC_MALLOC(sizeof(PID_TYPE)*dfc->CompactTable4[i].array[n].cnt,
																							DFC_MEMORY_TYPE__CT4);	
				memcpy(tempPID, dfc->CompactTable4[i].array[n].pid, 
														sizeof(PID_TYPE) * dfc->CompactTable4[i].array[n].cnt);
				//free(dfc->CompactTable4[i].array[n].pid);
				DFC_FREE(dfc->CompactTable4[i].array[n].pid, dfc->CompactTable4[i].array[n].cnt*sizeof(PID_TYPE),
															DFC_MEMORY_TYPE__CT4);
				dfc->CompactTable4[i].array[n].pid = NULL;

				int temp_cnt = 0; // cnt for 4 byte patterns.

				for(m = 0; m < dfc->CompactTable4[i].array[n].cnt; m++){
					//DFC_PATTERN *mlist = dfc->dfcMatchList[]; 
					int pat_len = dfc->dfcMatchList[tempPID[m]]->n - 4;

					if(pat_len == 0){ /* When pat length is 4 */
						temp_cnt ++;
						dfc->CompactTable4[i].array[n].pid =
							 (PID_TYPE *)realloc(dfc->CompactTable4[i].array[n].pid, sizeof(PID_TYPE)*temp_cnt);
						dfc->CompactTable4[i].array[n].pid[temp_cnt-1] = tempPID[m];
					}else if (pat_len == 1){ /* When pat length is 5 */
						if(dfc->dfcMatchList[tempPID[m]]->nocase){
							temp[1] = tolower(dfc->dfcMatchList[tempPID[m]]->patrn[0]);
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable4[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable4[i].array[n].CompactTable, temp,	tempPID[m],
																							DFC_MEMORY_TYPE__CT4);
							}

							temp[1] = toupper(dfc->dfcMatchList[tempPID[m]]->patrn[0]);
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable4[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable4[i].array[n].CompactTable, temp,	tempPID[m],
																							DFC_MEMORY_TYPE__CT4);
							}
						}else{
							temp[1] = dfc->dfcMatchList[tempPID[m]]->casepatrn[0];
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable4[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable4[i].array[n].CompactTable, temp,	tempPID[m],
																						DFC_MEMORY_TYPE__CT4);
							}
						}
						
					}else { /* When pat length is 7 (pat_len is equal to 2(6) or 3(7)) */
						if(dfc->dfcMatchList[tempPID[m]]->nocase){
							alpha_cnt = 0;
							do{
								for (l=1, k = 0; l >=0; --l, k++) {
									flag[k] =  (alpha_cnt >> l) & 1;
								}

								for (l = pat_len - 2, k = 0; l <= pat_len - 1; l++, k++){
									Build_pattern(dfc->dfcMatchList[tempPID[m]], flag, temp, 0,l,k);
								}

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable4[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable4[i].array[n].CompactTable,temp,tempPID[m],
																						DFC_MEMORY_TYPE__CT4);

								alpha_cnt++;
							}while(alpha_cnt < 4);
						}else{ /* case sensitive pattern */
							temp[0] = dfc->dfcMatchList[tempPID[m]]->casepatrn[pat_len - 2];
							temp[1] = dfc->dfcMatchList[tempPID[m]]->casepatrn[pat_len - 1];

                            fragment_16 = (temp[1] << 8) | temp[0];
                            byteIndex = BINDEX(fragment_16);
                            bitMask = BMASK(fragment_16);

							dfc->CompactTable4[i].array[n].DirectFilter[byteIndex] |= bitMask;

							Add_PID_to_2B_CT(dfc->CompactTable4[i].array[n].CompactTable,temp,tempPID[m],
																						 DFC_MEMORY_TYPE__CT4);
						}
					}
				}

				dfc->CompactTable4[i].array[n].cnt = temp_cnt;
				DFC_FREE(tempPID, sizeof(PID_TYPE)*dfc->CompactTable4[i].array[n].cnt, DFC_MEMORY_TYPE__CT4);
			}
		}
	}

#if 1
	/* For CT8 */
	for(i=0; i < CT8_TABLE_SIZE; i++){
		for(n=0; n < dfc->CompactTable8[i].cnt; n++){
			/* If the number of PID is bigger than RECURSIVE_BOUNDARY, do recursive filtering */
			if(dfc->CompactTable8[i].array[n].cnt >= RECURSIVE_BOUNDARY){
				/* Initialization */
				dfc->CompactTable8[i].array[n].DirectFilter
                    = (uint8_t*)DFC_MALLOC(DF_SIZE_REAL*sizeof(uint8_t), DFC_MEMORY_TYPE__CT8);
				dfc->CompactTable8[i].array[n].CompactTable=(CT_Type_2_2B *)DFC_MALLOC
                                                (sizeof(CT_Type_2_2B)*RECURSIVE_CT_SIZE,DFC_MEMORY_TYPE__CT8);

				if(!dfc->CompactTable8[i].array[n].DirectFilter || !dfc->CompactTable8[i].array[n].CompactTable){
					printe("Failed to allocate memory for recursive things.\n");
					return -1;
				}

				PID_TYPE *tempPID = (PID_TYPE*)DFC_MALLOC(sizeof(PID_TYPE)*dfc->CompactTable8[i].array[n].cnt,
																						DFC_MEMORY_TYPE__CT8);	
				memcpy(tempPID, dfc->CompactTable8[i].array[n].pid, 
														sizeof(PID_TYPE) * dfc->CompactTable8[i].array[n].cnt);
				//free(dfc->CompactTable8[i].array[n].pid);
				DFC_FREE(dfc->CompactTable8[i].array[n].pid, dfc->CompactTable8[i].array[n].cnt*sizeof(PID_TYPE),
															DFC_MEMORY_TYPE__CT8);
				dfc->CompactTable8[i].array[n].pid = NULL;

				int temp_cnt = 0; // cnt for 8 byte patterns.

				for(m = 0; m < dfc->CompactTable8[i].array[n].cnt; m++){
					//DFC_PATTERN *mlist = dfc->dfcMatchList[]; 
					int pat_len = dfc->dfcMatchList[tempPID[m]]->n - 8;

					if(pat_len == 0){ /* When pat length is 8 */
						temp_cnt ++;
						dfc->CompactTable8[i].array[n].pid =
							 (PID_TYPE *)realloc(dfc->CompactTable8[i].array[n].pid, sizeof(PID_TYPE)*temp_cnt);
						dfc->CompactTable8[i].array[n].pid[temp_cnt-1] = tempPID[m];
					}else if (pat_len == 1){ /* When pat length is 9 */
						if(dfc->dfcMatchList[tempPID[m]]->nocase){
							temp[1] = tolower(dfc->dfcMatchList[tempPID[m]]->patrn[0]);
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable8[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable8[i].array[n].CompactTable, temp,	tempPID[m],
																						DFC_MEMORY_TYPE__CT8);
							}

							temp[1] = toupper(dfc->dfcMatchList[tempPID[m]]->patrn[0]);
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable8[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable8[i].array[n].CompactTable, temp,	tempPID[m],
																							DFC_MEMORY_TYPE__CT8);
							}
						}else{
							temp[1] = dfc->dfcMatchList[tempPID[m]]->casepatrn[0];
							for(l = 0; l < 256; l++){
								temp[0] = l;

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable8[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable8[i].array[n].CompactTable, temp,	tempPID[m],
																						DFC_MEMORY_TYPE__CT8);
							}
						}
						
					}else { /* longer than or equal to 10 */
						if(dfc->dfcMatchList[tempPID[m]]->nocase){
							alpha_cnt = 0;
							do{
								for (l=1, k = 0; l >=0; --l, k++) {
									flag[k] =  (alpha_cnt >> l) & 1;
								}

								for (l = pat_len - 2, k = 0; l <= pat_len - 1; l++, k++){
									Build_pattern(dfc->dfcMatchList[tempPID[m]], flag, temp, 0,l,k);
								}

                                fragment_16 = (temp[1] << 8) | temp[0];
                                byteIndex = BINDEX(fragment_16);
                                bitMask = BMASK(fragment_16);

								dfc->CompactTable8[i].array[n].DirectFilter[byteIndex] |= bitMask;

								Add_PID_to_2B_CT(dfc->CompactTable8[i].array[n].CompactTable,temp,tempPID[m],
																						DFC_MEMORY_TYPE__CT8);

								alpha_cnt++;
							}while(alpha_cnt < 4);
						}else{ /* case sensitive pattern */
							temp[0] = dfc->dfcMatchList[tempPID[m]]->casepatrn[pat_len - 2];
							temp[1] = dfc->dfcMatchList[tempPID[m]]->casepatrn[pat_len - 1];

                            fragment_16 = (temp[1] << 8) | temp[0];
                            byteIndex = BINDEX(fragment_16);
                            bitMask = BMASK(fragment_16);

							dfc->CompactTable8[i].array[n].DirectFilter[byteIndex] |= bitMask;

							Add_PID_to_2B_CT(dfc->CompactTable8[i].array[n].CompactTable,temp,tempPID[m],
																							DFC_MEMORY_TYPE__CT8);
						}
					}
				}

				dfc->CompactTable8[i].array[n].cnt = temp_cnt;
				DFC_FREE(tempPID, sizeof(PID_TYPE)*dfc->CompactTable8[i].array[n].cnt, DFC_MEMORY_TYPE__CT8);
			}
		}
	}
#endif
#endif

/* ####################################################################################### */

    //printf("Recursive Initialization is done.\n");

/* ####################################################################################### */
/* ###############                   Print Information                    ################ */
/* ####################################################################################### */
#ifdef PRINT_INFO
	//DFC_PrintInfo(dfc);
#endif
/* ####################################################################################### */

	return 0;
}

static int Verification_CT1(VERIFI_ARGUMENT)
{
	int i;
	for(i = 0; i < dfc->CompactTable1[*(buf-2)].cnt; i++){
		PID_TYPE pid = dfc->CompactTable1[*(buf-2)].pid[i];
		DFC_PATTERN *mlist = dfc->dfcMatchList[pid];

        ACTION_FOR_MATCH;
	}
	return matches;
}

static int Verification_CT2(VERIFI_ARGUMENT)
{
	uint32_t crc = crc64(0, (uint64_t)(*(uint16_t*)(buf - 2)), 2); //_mm_crc32_u16(0, *(uint16_t*)(buf-2));

	// 2. calculate index
	crc &= CT2_TABLE_SIZE_MASK;

	BUC_CNT_TYPE i;
	for(i = 0; i < dfc->CompactTable2[crc].cnt; i++){
		if(dfc->CompactTable2[crc].array[i].pat == *(uint16_t*)(buf-2)){
			PID_CNT_TYPE j;
			if(!dfc->CompactTable2[crc].array[i].DirectFilter){
				for(j = 0; j < dfc->CompactTable2[crc].array[i].cnt; j++){
					PID_TYPE pid = dfc->CompactTable2[crc].array[i].pid[j];
					
					DFC_PATTERN *mlist = dfc->dfcMatchList[pid];
					if(buf - starting_point >= mlist->n){
						if(mlist->nocase){
							if(!my_strncasecmp(buf-(mlist->n), mlist->casepatrn, mlist->n - 2)){
                                ACTION_FOR_MATCH;
							}
						}else{
							if(!my_strncmp(buf-(mlist->n), mlist->casepatrn, mlist->n - 2)){
                                ACTION_FOR_MATCH;
							}
						}
					}
				}
			}else{
				for(j = 0; j < dfc->CompactTable2[crc].array[i].cnt; j++){
                    PID_TYPE pid = dfc->CompactTable2[crc].array[i].pid[j];
                    DFC_PATTERN *mlist = dfc->dfcMatchList[pid];

                    ACTION_FOR_MATCH;
				}
				DTYPE data = *(uint16_t*)(buf-4);
				BTYPE index = BINDEX(data);
				BTYPE mask = BMASK(data);	

				if(dfc->CompactTable2[crc].array[i].DirectFilter[index] & mask) {
					uint32_t crc2 = crc64(0, (uint64_t)data, 2); //_mm_crc32_u16(0, data);

					// 2. calculate index
					crc2 &= CT2_TABLE_SIZE_MASK;

					BUC_CNT_TYPE k;
					for(k = 0; k < dfc->CompactTable2[crc].array[i].CompactTable[crc2].cnt; k++){
						if(dfc->CompactTable2[crc].array[i].CompactTable[crc2].array[k].pat == data){
							PID_CNT_TYPE l;
							for(l = 0; l < dfc->CompactTable2[crc].array[i].CompactTable[crc2].array[k].cnt; l++){
								PID_TYPE pid = dfc->CompactTable2[crc].array[i].CompactTable[crc2].array[k].pid[l];
								DFC_PATTERN *mlist = dfc->dfcMatchList[pid];

                                ACTION_FOR_MATCH;
							}
							break;
						}
					}
				}
			}
			break;
		}
	}
	return matches;
}


static int Verification_CT4_7(VERIFI_ARGUMENT)
{
	// 1. Convert payload to uppercase
	unsigned char *temp = buf-2;

	// 2. calculate crc
	uint32_t crc = crc64(0, (uint64_t)(*(uint32_t*)temp), 4); //_mm_crc32_u32(0, *(uint32_t*)temp);

	// 3. calculate index
	crc &= CT4_TABLE_SIZE_MASK;

	// 4.
	BUC_CNT_TYPE i;
	for(i = 0; i < dfc->CompactTable4[crc].cnt; i++){
		if(dfc->CompactTable4[crc].array[i].pat == *(uint32_t*)temp){
			PID_CNT_TYPE j;
			if(!dfc->CompactTable4[crc].array[i].DirectFilter){
				for(j = 0; j < dfc->CompactTable4[crc].array[i].cnt; j++){
					PID_TYPE pid = dfc->CompactTable4[crc].array[i].pid[j];
					
					DFC_PATTERN *mlist = dfc->dfcMatchList[pid];
					if(buf - starting_point >= mlist->n - 2){
						if(mlist->nocase){
							if(!my_strncasecmp(buf-(mlist->n-2), mlist->casepatrn, mlist->n - 4)){
                                ACTION_FOR_MATCH;
							}
						}else{
							if(!my_strncmp(buf-(mlist->n-2), mlist->casepatrn, mlist->n - 4)){
                                ACTION_FOR_MATCH;
							}
						}
					}
				}
			}else{
				for(j = 0; j < dfc->CompactTable4[crc].array[i].cnt; j++){
                    PID_TYPE pid = dfc->CompactTable4[crc].array[i].pid[j];
                    DFC_PATTERN *mlist = dfc->dfcMatchList[pid];

                    ACTION_FOR_MATCH;
				}

				DTYPE data = *(uint16_t*)(buf-4);
				BTYPE index = BINDEX(data);
				BTYPE mask = BMASK(data);	

				if(dfc->CompactTable4[crc].array[i].DirectFilter[index] & mask) {
					uint32_t crc2 = crc64(0, (uint64_t)data, 2); //_mm_crc32_u16(0, data);

					// 2. calculate index
					crc2 &= CT2_TABLE_SIZE_MASK;

					BUC_CNT_TYPE k;
					for(k = 0; k < dfc->CompactTable4[crc].array[i].CompactTable[crc2].cnt; k++){
						if(dfc->CompactTable4[crc].array[i].CompactTable[crc2].array[k].pat == data){
							PID_CNT_TYPE l;
							for(l = 0; l < dfc->CompactTable4[crc].array[i].CompactTable[crc2].array[k].cnt; l++){
								PID_TYPE pid = dfc->CompactTable4[crc].array[i].CompactTable[crc2].array[k].pid[l];
								
								DFC_PATTERN *mlist = dfc->dfcMatchList[pid];
								if(mlist->nocase){
									if(!my_strncasecmp(buf-(mlist->n-2), mlist->casepatrn, mlist->n - 6)){
                                        ACTION_FOR_MATCH;
									}
								}else{
									if(!my_strncmp(buf-(mlist->n-2), mlist->casepatrn, mlist->n - 6)){
                                        ACTION_FOR_MATCH;
									}
								}
							}
							break;
						}
					}
				}
			}
			break;
		}
	}
	return matches;
}



static int Verification_CT8_plus(VERIFI_ARGUMENT)
{
	// 1. Convert payload to uppercase
#ifdef CT8_SWITCH
	unsigned char temp[8];
	ConvertCaseEx (temp, buf-2, 8);
#else
	unsigned char *temp = buf-2;
#endif

	// 2. calculate crc
    uint32_t fragment_32 = (temp[7] << 24) | (temp[6] << 16) | (temp[5] << 8) | temp[4];
    uint64_t fragment_64 = ((uint64_t)fragment_32 << 32) |
                           (temp[3] << 24) | (temp[2] << 16) | (temp[1] << 8) | temp[0];
	uint64_t crc = crc64(0, (uint64_t)fragment_64, 8); //_mm_crc32_u64(0, fragment_64);

	// 3. calculate index
	crc &= CT8_TABLE_SIZE_MASK;

	BUC_CNT_TYPE i;
	for(i = 0; i < dfc->CompactTable8[crc].cnt; i++){
		if(dfc->CompactTable8[crc].array[i].pat == fragment_64){
			//matches++;	break;

			PID_CNT_TYPE j;
			if(!dfc->CompactTable8[crc].array[i].DirectFilter){
				for(j = 0; j < dfc->CompactTable8[crc].array[i].cnt; j++){
					PID_TYPE pid = dfc->CompactTable8[crc].array[i].pid[j];
					
					DFC_PATTERN *mlist = dfc->dfcMatchList[pid];
					int comparison_requirement = min_pattern_interval*(mlist->n-8)/pattern_interval+2;
					if(buf - starting_point >= comparison_requirement){
						if(mlist->nocase){
							if(!my_strncasecmp(buf-comparison_requirement, mlist->casepatrn, mlist->n)){
                                ACTION_FOR_MATCH;
							}
						}else{
							if(!my_strncmp(buf-comparison_requirement, mlist->casepatrn, mlist->n)){
                                ACTION_FOR_MATCH;
							}
						}
					}
				}
#if 1
			}else{
				for(j = 0; j < dfc->CompactTable8[crc].array[i].cnt; j++){
                    PID_TYPE pid = dfc->CompactTable8[crc].array[i].pid[j];

                    DFC_PATTERN *mlist = dfc->dfcMatchList[pid];
					int comparison_requirement = min_pattern_interval*(mlist->n-8)/pattern_interval+2;
					if(mlist->nocase){
						if(!my_strncasecmp(buf-comparison_requirement, mlist->casepatrn, mlist->n)){
                            ACTION_FOR_MATCH;
						}
					}else{
						if(!my_strncmp(buf-comparison_requirement, mlist->casepatrn, mlist->n)){
                            ACTION_FOR_MATCH;
						}
					}
				}

				DTYPE data = *(uint16_t*)(buf-4);
				BTYPE index = BINDEX(data);
				BTYPE mask = BMASK(data);	

				if(dfc->CompactTable8[crc].array[i].DirectFilter[index] & mask) {
					uint32_t crc2 = crc64(0, (uint64_t)data, 2); //_mm_crc32_u16(0, data);

					// 2. calculate index
					crc2 &= CT2_TABLE_SIZE_MASK;

					BUC_CNT_TYPE k;
					for(k = 0; k < dfc->CompactTable8[crc].array[i].CompactTable[crc2].cnt; k++){
						if(dfc->CompactTable8[crc].array[i].CompactTable[crc2].array[k].pat == data){
							PID_CNT_TYPE l;
							for(l = 0; l < dfc->CompactTable8[crc].array[i].CompactTable[crc2].array[k].cnt; l++){
								PID_TYPE pid = dfc->CompactTable8[crc].array[i].CompactTable[crc2].array[k].pid[l];
								
								DFC_PATTERN *mlist = dfc->dfcMatchList[pid];
								int comparison_requirement = min_pattern_interval*(mlist->n-8)/pattern_interval+2;
								if(buf - starting_point >= comparison_requirement){
									if(mlist->nocase){
										if(!my_strncasecmp(buf-comparison_requirement, mlist->casepatrn,
																					 mlist->n)){
                                            ACTION_FOR_MATCH;
										}
									}else{
										if(!my_strncmp(buf-comparison_requirement, mlist->casepatrn,
																					 mlist->n)){
                                            ACTION_FOR_MATCH;
										}
									}

								}
							}
							break;
						}
					}
				}
			}
#endif
			break;
		}
	}

	return matches;
}


static inline int Progressive_Filtering(PROGRE_ARGUMENT) 
{
#ifdef ENABLE_PROGRESSIVE_FILTERING
    if(dfc->cDF0[*(buf-2)]){
		matches = Verification_CT1(VERIFI_PARAMETER);
	}

    if (unlikely(dfc->cDF1[idx] & msk)) {
		matches = Verification_CT2(VERIFI_PARAMETER);
	}

    if(rest_len >= 4){
        DTYPE data = *(uint16_t*)(buf);
        BTYPE index = BINDEX(data);
        BTYPE mask = BMASK(data);	

        if (unlikely((mask & dfc->ADD_DF_4_plus[index]))) {
            if (unlikely(mask & dfc->ADD_DF_4_1[index])) {
                //if (unlikely(msk & dfc->cDF2[idx])) {
                    matches = Verification_CT4_7(VERIFI_PARAMETER);
                    //matches ++;
                //}
            }

            DTYPE data8 = *(uint16_t*)(&buf[4]);
            BTYPE index8 = BINDEX(data8);
            BTYPE mask8 = BMASK(data8);

            if (unlikely(mask8 & dfc->ADD_DF_8_1[index8])) {
                data8 = *(uint16_t*)(&buf[2]);
                index8 = BINDEX(data8);
                mask8 = BMASK(data8);
                if (unlikely(mask8 & dfc->ADD_DF_8_2[index8])) {
                    if ((rest_len >= 8)) {
                        matches = Verification_CT8_plus(VERIFI_PARAMETER);
                        //matches ++;
                    }
                }
            }
        }
    }
#else
	matches = Verification_CT1(VERIFI_PARAMETER);
	matches = Verification_CT2(VERIFI_PARAMETER);
	matches = Verification_CT4_7(VERIFI_PARAMETER);
	matches = Verification_CT8_plus(VERIFI_PARAMETER);
#endif

    return matches;
}

int DFC_Search(SEARCH_ARGUMENT)
{
	int i;
    int matches = 0;

	if (unlikely(buflen <= 0))
        return 0;

	uint8_t *DirectFilter1 = dfc->DirectFilter1;

	for (i = 0; i < buflen - 1; i++) {
		DTYPE data = *(uint16_t*)(&buf[i]);
		BTYPE index = BINDEX(data);
		BTYPE mask = BMASK(data);

		if (unlikely(DirectFilter1[index] & mask)) {
			matches = Progressive_Filtering(PROGRE_PARAMETER);
		}
	} 

	/* It is needed to check last 1 byte from payload */
    if(dfc->cDF0[buf[buflen-1]]){
        int i;
        for(i = 0; i < dfc->CompactTable1[buf[buflen-1]].cnt; i++){
            PID_TYPE pid = dfc->CompactTable1[buf[buflen-1]].pid[i];
			DFC_PATTERN *mlist = dfc->dfcMatchList[pid];

            ACTION_FOR_MATCH;
        }
    }

	return matches;
}


/*************************************************************************************/
/*                                       Utility                                     */
/*************************************************************************************/
static float my_sqrtf(float input, float x)
{
    int i;
	if( x == 0 && input == 0)
		return 0;
    for(i = 0; i < 10; i++){
        x = (x + (input / x))/2;
    }
    return x;
}

static void init_xlatcase()
{
	int i;
	for (i = 0; i < 256; i++)
		xlatcase[i] = (unsigned char)toupper(i);
}

static inline void ConvertCaseEx (unsigned char *d, unsigned char *s, int m)
{
	int i;

	for (i=0; i < m; i++)
		d[i] = xlatcase[ s[i] ];
}

static inline int my_strncmp(unsigned char *a, unsigned char *b, int n){
    int i;
    for(i = 0; i < n ; i++){
        if(a[i] != b[i])
            return -1;
    }
    return 0;
}

static inline int my_strncasecmp(unsigned char *a, unsigned char *b, int n){
    int i;
    for(i = 0; i < n ; i++){
        if(tolower(a[i]) != tolower(b[i]))
            return -1;
    }
    return 0;
}

static void * DFC_REALLOC(void *p, uint16_t n, dfcDataType type, dfcMemoryType type2 )
{
	switch(type){
		case DFC_PID_TYPE:
			p = realloc((PID_TYPE*)p, sizeof(PID_TYPE)*n);
			dfc_total_memory += sizeof(PID_TYPE);
			switch(type2){
                case DFC_MEMORY_TYPE__PATTERN:
                    dfc_pattern_memory += sizeof(PID_TYPE);
                    break;
				case DFC_MEMORY_TYPE__CT2:
					dfc_memory_ct2 += sizeof(PID_TYPE);
					break;
				case DFC_MEMORY_TYPE__CT3:
					dfc_memory_ct3 += sizeof(PID_TYPE);
					break;
				case DFC_MEMORY_TYPE__CT4:
					dfc_memory_ct4 += sizeof(PID_TYPE);
					break;
				case DFC_MEMORY_TYPE__CT8:
					dfc_memory_ct8 += sizeof(PID_TYPE);
					break;
				default:
					break;
			}
			return p;
		case DFC_CT_Type_2_Array:
			p = realloc((CT_Type_2_Array*)p, sizeof(CT_Type_2_Array)*n);
			dfc_total_memory += sizeof(CT_Type_2_Array);
			switch(type2){
				case DFC_MEMORY_TYPE__CT2:
					dfc_memory_ct2 += sizeof(CT_Type_2_Array);
					break;
				case DFC_MEMORY_TYPE__CT3:
					dfc_memory_ct3 += sizeof(CT_Type_2_Array);
					break;
				case DFC_MEMORY_TYPE__CT4:
					dfc_memory_ct4 += sizeof(CT_Type_2_Array);
					break;
				case DFC_MEMORY_TYPE__CT8:
					dfc_memory_ct8 += sizeof(CT_Type_2_Array);
					break;
				default:
					break;
			}
			return p;
		case DFC_CT_Type_2_2B_Array:
			p = realloc((CT_Type_2_2B_Array*)p, sizeof(CT_Type_2_2B_Array)*n);
			dfc_total_memory += sizeof(CT_Type_2_2B_Array);
			switch(type2){
				case DFC_MEMORY_TYPE__CT2:
					dfc_memory_ct2 += sizeof(CT_Type_2_2B_Array);
					break;
				case DFC_MEMORY_TYPE__CT3:
					dfc_memory_ct3 += sizeof(CT_Type_2_2B_Array);
					break;
				case DFC_MEMORY_TYPE__CT4:
					dfc_memory_ct4 += sizeof(CT_Type_2_2B_Array);
					break;
				case DFC_MEMORY_TYPE__CT8:
					dfc_memory_ct8 += sizeof(CT_Type_2_2B_Array);
					break;
				default:
					break;
			}
			return p;
		case DFC_CT_Type_2_8B_Array:
			p = realloc((CT_Type_2_8B_Array*)p, sizeof(CT_Type_2_8B_Array)*n);
			dfc_total_memory += sizeof(CT_Type_2_8B_Array);
			switch(type2){
				case DFC_MEMORY_TYPE__CT2:
					dfc_memory_ct2 += sizeof(CT_Type_2_8B_Array);
					break;
				case DFC_MEMORY_TYPE__CT3:
					dfc_memory_ct3 += sizeof(CT_Type_2_8B_Array);
					break;
				case DFC_MEMORY_TYPE__CT4:
					dfc_memory_ct4 += sizeof(CT_Type_2_8B_Array);
					break;
				case DFC_MEMORY_TYPE__CT8:
					dfc_memory_ct8 += sizeof(CT_Type_2_8B_Array);
					break;
				default:
					break;
			}
			return p;
		default:
			printf("ERROR! Data Type is not correct!\n");
			break;
	}
	return NULL;
}

static void DFC_FREE(void *p, int n, dfcMemoryType type){
	free(p);
	switch (type)
	{
		case DFC_MEMORY_TYPE__DFC:
			break;
		case DFC_MEMORY_TYPE__PATTERN:
			dfc_pattern_memory -= n;
			break;
		case DFC_MEMORY_TYPE__CT2:
			dfc_memory_ct2 -= n;
			break;
		case DFC_MEMORY_TYPE__CT3:
			dfc_memory_ct3 -= n;
			break;
		case DFC_MEMORY_TYPE__CT4:
			dfc_memory_ct4 -= n;
			break;
		case DFC_MEMORY_TYPE__CT8:
			dfc_memory_ct8 -= n;
			break;
		case DFC_MEMORY_TYPE__NONE:
			break;
		default:
			//printf("%s(%d) Invalid memory type\n", __FILE__, __LINE__);
			break;
	}
	dfc_total_memory -= n;
}

static void * DFC_MALLOC(int n, dfcMemoryType type )
{
    void *p = calloc(1, n); // initialize it to 0

    if (p != NULL)
    {
        switch (type)
        {
            case DFC_MEMORY_TYPE__DFC:
                break;
            case DFC_MEMORY_TYPE__PATTERN:
				dfc_pattern_memory += n;
                break;
            case DFC_MEMORY_TYPE__CT2:
				dfc_memory_ct2 += n;
                break;
            case DFC_MEMORY_TYPE__CT3:
				dfc_memory_ct3 += n;
                break;
            case DFC_MEMORY_TYPE__CT4:
				dfc_memory_ct4 += n;
                break;
            case DFC_MEMORY_TYPE__CT8:
				dfc_memory_ct8 += n;
                break;
            case DFC_MEMORY_TYPE__NONE:
                break;
            default:
                printf("%s(%d) Invalid memory type\n", __FILE__, __LINE__);
                break;
        }
        dfc_total_memory += n;
    }
    return p;
}

static void Build_pattern(DFC_PATTERN *p, uint8_t *flag, uint8_t *temp, uint32_t i, int j, int k){
    if(p->nocase){
        if((p->patrn[j] >= 65 && p->patrn[j]<=90) ||
            (p->patrn[j] >= 97 && p->patrn[j]<=122)) {
            if(flag[k] == 0)
                temp[k] = tolower(p->patrn[j]);
            else
                temp[k] = toupper(p->patrn[j]);
        } else {
            temp[k] = p->patrn[j];
        }
    }else{
        temp[k] = p->casepatrn[j]; // original pattern
    }
    return ;
}

/**
 * \internal
 * \brief Creates a hash of the pattern.  We use it for the hashing process
 *        during the initial pattern insertion time, to cull duplicate sigs.
 *
 * \param pat    Pointer to the pattern.
 * \param patlen Pattern length.
 *
 * \retval hash A 32 bit unsigned hash.
 */
static inline uint32_t DFC_InitHashRaw(uint8_t *pat, uint16_t patlen, int nocase)
{
    uint32_t hash  = 0;
    int i;
    for (i = 0; i < patlen; i++)
        hash+= pat[i];

    hash += nocase;

    return (hash % INIT_HASH_SIZE);
}


/**
 * \internal
 * \brief Looks up a pattern.  We use it for the hashing process during the
 *        the initial pattern insertion time, to cull duplicate sigs.
 *
 * \param ctx    Pointer to the DFC structure.
 * \param pat    Pointer to the pattern.
 * \param patlen Pattern length.
 * \param nocase Pattern case insensitivity
 *
 */
static inline DFC_PATTERN *DFC_InitHashLookup(DFC_STRUCTURE *ctx, uint8_t *pat,
                                              uint16_t patlen, int nocase)
{
    uint32_t hash = DFC_InitHashRaw(pat, patlen, nocase);

    if (ctx->init_hash == NULL) {
        return NULL;
    }

    DFC_PATTERN *t = ctx->init_hash[hash];
    for ( ; t != NULL; t = t->next) {
        if (!strcmp((char*)t->casepatrn, (char*)pat))
            return t;
    }

    return NULL;
}

static inline int DFC_InitHashAdd(DFC_STRUCTURE *ctx, DFC_PATTERN *p)
{
    uint32_t hash = DFC_InitHashRaw(p->casepatrn, p->n, p->nocase);

    if (ctx->init_hash == NULL) {
        return 0;
    }

    if (ctx->init_hash[hash] == NULL) {
        ctx->init_hash[hash] = p;
        return 0;
    }

    DFC_PATTERN *tt = NULL;
    DFC_PATTERN *t = ctx->init_hash[hash];

    /* get the list tail */
    do {
        tt = t;
        t = t->next;
    } while (t != NULL);

    tt->next = p;

    return 0;
}

/*************************************************************************************/
//}
//#endif
