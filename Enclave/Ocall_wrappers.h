#ifndef _OCALL_WRAPPERS_H_
#define _OCALL_WRAPPERS_H_

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

#include "sgx_trts.h"

#include "Enclave_t.h"

#if defined(__cplusplus)
extern "C" {
#endif

int sgx_printf(const char *fmt, ...);
void sgx_printe(const char *fname, const char *fmt, ...);
void sgx_printl(const char *fname, const char *fmt, ...);
void sgx_exit(int exit_status);

#define printf sgx_printf
#define printe(fmt, ...) sgx_printe(__FUNCTION__, fmt, ##__VA_ARGS__)
#define printl(fmt, ...) sgx_printl(__FUNCTION__, fmt, ##__VA_ARGS__)

#if defined(__cplusplus)
}
#endif

#endif /* !_OCALL_WRAPPERS_H_ */