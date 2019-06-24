#ifndef _OCALL_IMPLEMENTS_H_
#define _OCALL_IMPLEMENTS_H_

#include <stdio.h>
#include <string.h>

#if defined(__cplusplus)
extern "C" {
#endif

void ocall_print_string(const char *str);
void ocall_sgx_exit(int e);

#if defined(__cplusplus)
}
#endif

#endif /* !_OCALL_IMPLEMENTS_H_ */