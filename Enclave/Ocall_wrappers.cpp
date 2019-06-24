#include "Ocall_wrappers.h"

int sgx_printf(const char *fmt, ...)
{
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    ocall_print_string(buf);
}

void sgx_printe(const char *fname, const char *fmt, ...)
{
    char ebuf[BUFSIZ] = {'\0'};
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    snprintf(ebuf, sizeof(ebuf), "Error: %s failed!: %s\n", fname, buf);
    ocall_print_string(ebuf);
}

void sgx_printl(const char *fname, const char *fmt, ...)
{
    char ebuf[BUFSIZ] = {'\0'};
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    snprintf(ebuf, sizeof(ebuf), "LOG: %s: %s\n", fname, buf);
    ocall_print_string(ebuf);
}

void sgx_exit(int exit_status)
{
    printf("sgx_exit: exit(%d) called!\n", exit_status);
    ocall_sgx_exit(exit_status);
}