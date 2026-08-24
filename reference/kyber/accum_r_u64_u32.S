    .text
    .globl brlt_accum_r_u64_u32_asm
    .type brlt_accum_r_u64_u32_asm, @function

// void brlt_accum_r_u64_u32_asm(uint64_t *acc, const uint32_t *src, uint32_t n, uint64_t r)
// SysV x86_64: rdi, rsi, edx, rcx
brlt_accum_r_u64_u32_asm:
    testl   %edx, %edx
    je      .Ldone

.Lloop4:
    cmpl    $4, %edx
    jb      .Ltail

    movl    0(%rsi), %r8d
    imulq   %rcx, %r8
    addq    %r8, 0(%rdi)

    movl    4(%rsi), %r9d
    imulq   %rcx, %r9
    addq    %r9, 8(%rdi)

    movl    8(%rsi), %r10d
    imulq   %rcx, %r10
    addq    %r10, 16(%rdi)

    movl    12(%rsi), %r11d
    imulq   %rcx, %r11
    addq    %r11, 24(%rdi)

    addq    $16, %rsi
    addq    $32, %rdi
    subl    $4, %edx
    jmp     .Lloop4

.Ltail:
    testl   %edx, %edx
    je      .Ldone

.Lloop1:
    movl    (%rsi), %r8d
    imulq   %rcx, %r8
    addq    %r8, (%rdi)

    addq    $4, %rsi
    addq    $8, %rdi
    decl    %edx
    jne     .Lloop1

.Ldone:
    ret

    .size brlt_accum_r_u64_u32_asm, .-brlt_accum_r_u64_u32_asm
