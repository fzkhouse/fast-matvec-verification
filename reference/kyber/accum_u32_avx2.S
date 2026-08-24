    .text
    .globl brlt_accum_u32_avx2
    .type brlt_accum_u32_avx2, @function

// void brlt_accum_u32_avx2(uint32_t *acc, const uint32_t *src, uint32_t n)
brlt_accum_u32_avx2:
    testl   %edx, %edx
    je      .Ldone

    movl    %edx, %ecx

.Lloop16:
    cmpl    $16, %ecx
    jb      .Lloop8

    vmovdqu (%rdi), %ymm0
    vmovdqu 32(%rdi), %ymm1
    vmovdqu (%rsi), %ymm2
    vmovdqu 32(%rsi), %ymm3

    vpaddd  %ymm2, %ymm0, %ymm0
    vpaddd  %ymm3, %ymm1, %ymm1

    vmovdqu %ymm0, (%rdi)
    vmovdqu %ymm1, 32(%rdi)

    addq    $64, %rdi
    addq    $64, %rsi
    subl    $16, %ecx
    jmp     .Lloop16

.Lloop8:
    cmpl    $8, %ecx
    jb      .Ltail

    vmovdqu (%rdi), %ymm0
    vmovdqu (%rsi), %ymm1
    vpaddd  %ymm1, %ymm0, %ymm0
    vmovdqu %ymm0, (%rdi)

    addq    $32, %rdi
    addq    $32, %rsi
    subl    $8, %ecx

.Ltail:
    testl   %ecx, %ecx
    je      .Ldone_scalar

.Lscalar:
    movl    (%rsi), %eax
    addl    %eax, (%rdi)
    addq    $4, %rdi
    addq    $4, %rsi
    decl    %ecx
    jne     .Lscalar

.Ldone_scalar:
    vzeroupper

.Ldone:
    ret

    .size brlt_accum_u32_avx2, .-brlt_accum_u32_avx2
