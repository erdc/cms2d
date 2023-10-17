/*  ** CMS Preprocessor Definitions  **  */
/* Turn ON/OFF Unit Testing */
#undef UNIT_TEST

/* Turn ON/OFF Code Profiling */
#undef PROFILE

/* Turn ON/OFF Diagnostic Mode */
#undef DIAG_MODE

/* Turn ON/OFF Developmental Code */
#undef DEV_MODE

/* Checks */
#ifdef  _WIN32
#define   XMDF_IO 
#define   PROJ_CONV
#else
#undef    XMDF_IO
#undef    PROJ_CONV 
#endif
