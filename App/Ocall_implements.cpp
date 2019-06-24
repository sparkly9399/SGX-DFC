#include "Ocall_implements.h"
#include "Enclave_u.h"

void ocall_print_string(const char *str)
{
    printf("%s", str);
}

void ocall_sgx_exit(int e) 
{
    exit(e);
}