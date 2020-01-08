/*
** CMS Preprocessor Definitions
*/

/* Turn ON/OFF MergedCode flag */
#define MERGED_CODE

/* Turn ON/OFF Unit Testing */
#undef UNIT_TEST

/* Turn ON/OFF Code Profiling */
#undef PROFILE

/* Include XMDF Input and Output */
/* #define XMDF_IO  */

/* Turn ON/OFF Diagnostic Mode */
#undef DIAG_MODE

/* Turn ON/OFF Developmental Code */
#undef DEV_MODE

/* Turn ON/OFF Projection Conversion */
#define PROJ_CONV

/* Checks */
#ifdef  _WIN32
#define   XMDF_IO 
#else
#undef    XMDF_IO
#undef    PROJ_CONV 
#endif