function print_row(s1, s2, s3, s4, c)
    print("|", c)
    print(lpad(s1, 18, c))
    print(c, "|", c)
    print(lpad(s2, 10, c))
    print(c, "|", c)
    print(lpad(s3, 10, c))
    print(c, "|", c)
    print(lpad(s4, 10, c))
    return println(c, "|")
end
function prettyprint(name, ::Nothing)
    return print_row(name, "", "", "", ' ')
end
function prettyprint(name, b::BenchmarkTools.Trial)
    return print_row(
        name,
        BenchmarkTools.prettytime(time(b)),
        allocs(b),
        BenchmarkTools.prettymemory(memory(b)),
        ' ',
    )
end

function prettyprint(bs, bd, bt)
    print_row("", "Time", "Alloc", "Memory", ' ')
    print_row("", "", "", "", '-')
    prettyprint("SIMDPolynomials", bs)
    prettyprint("DynamicPolynomials", bd)
    return prettyprint("TypedPolynomials", bt)
end
