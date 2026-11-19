# Benchmark suite of FiniteFieldArithmetic.jl (BenchmarkTools `SUITE`, the convention of
# PkgBenchmark.jl and AirspeedVelocity.jl). Run it with `bench/run.jl`; results are recorded in
# `src/nimber.md`, "Measurements".
#
# - throughput: element-wise products of two arrays of N random elements (the bulk case, e.g.
#   learners), reported per element;
# - latency: a chain of N dependent products x ← x * y (each waits for the previous one);
# - baselines: XOR (field addition) and UInt8 integer multiplication over the same arrays.

using BenchmarkTools, Random
using FiniteFieldArithmetic
using FiniteFieldArithmetic: mul, inverse, Recursive, MulTable, LogTable, LogTable16, CLMul, bits

const N = 2^16              # elements per array for GF(256)
const N_LARGE = 2^12        # for GF(2^16) … GF(2^128), where Recursive is slow
const STRATEGIES = (Recursive(), MulTable(), LogTable(), LogTable16(), CLMul())
const LARGE = (("GF2_16", GF2_16), ("GF2_32", GF2_32), ("GF2_64", GF2_64), ("GF2_128", GF2_128))

# Elements per array (or chain steps) of each group, for reporting per element.
const GROUP_SIZE = Dict{String, Int}("GF256" => N, (name => N_LARGE for (name, _) in LARGE)...)

function chain(s, x, y, n)
    for _ in 1:n
        x = mul(s, x, y) + y        # "+ y" keeps x from collapsing to 0 or a short cycle
    end
    return x
end

const SUITE = BenchmarkGroup()

let rng = Xoshiro(2026)
    a, b = rand(rng, GF256, N), rand(rng, GF256, N)
    nz = [iszero(x) ? one(GF256) : x for x in b]
    c = similar(a)
    ua, ub, uc = bits.(a), bits.(b), similar(bits.(a))

    g = SUITE["GF256"] = BenchmarkGroup()
    for s in STRATEGIES
        name = string(nameof(typeof(s)))
        g["mul throughput", name] = @benchmarkable broadcast!((x, y) -> mul($s, x, y), $c, $a, $b)
        g["inv throughput", name] = @benchmarkable broadcast!(x -> inverse($s, x), $c, $nz)
        g["mul latency", name] = @benchmarkable chain($s, $(a[1]), $(b[1]), N)
    end
    g["mul throughput", "default (*)"] = @benchmarkable broadcast!(*, $c, $a, $b)
    g["baseline", "add (XOR)"] = @benchmarkable broadcast!(+, $c, $a, $b)
    g["baseline", "UInt8 *"] = @benchmarkable broadcast!(*, $uc, $ua, $ub)

    for (name, F) in LARGE
        a, b = rand(rng, F, N_LARGE), rand(rng, F, N_LARGE)
        nz = [iszero(x) ? one(F) : x for x in b]
        c = similar(a)
        ua, ub, uc = bits.(a), bits.(b), similar(bits.(a))
        g = SUITE[name] = BenchmarkGroup()
        for s in STRATEGIES
            sname = string(nameof(typeof(s)))
            g["mul throughput", sname] = @benchmarkable broadcast!((x, y) -> mul($s, x, y), $c, $a, $b)
            g["inv throughput", sname] = @benchmarkable broadcast!(x -> inverse($s, x), $c, $nz)
            g["mul latency", sname] = @benchmarkable chain($s, $(a[1]), $(b[1]), N_LARGE)
        end
        g["baseline", "add (XOR)"] = @benchmarkable broadcast!(+, $c, $a, $b)
        g["baseline", "integer *"] = @benchmarkable broadcast!(*, $uc, $ua, $ub)
    end
end

# --- Bulk and array operations (J8, src/bulk.md) ----------------------------------------------
#
# Per level, with N_BULK random elements: uniform sampling into an existing array, multiplication
# of an array by one constant (broadcast `c .* a` against the prepared `MulBy(c)`), and, for the
# levels a `BitslicedVector` supports, its element-wise product, packing and unpacking; all per
# element. Matrix products of two MATMUL_N × MATMUL_N matrices per strategy, per multiply-add
# (group "matmul <field>", MATMUL_N^3 multiply-adds).

using FiniteFieldArithmetic: ScalarMatMul, TableMatMul, CLMulMatMul, matmul!, mul_elementwise!

const N_BULK = 2^14
const MATMUL_N = 64
const BULK_LEVELS = (("GF256", GF256), LARGE...)
const BITSLICED_LEVELS = (GF256, GF2_16, GF2_32)
const MATMUL_STRATEGIES = (ScalarMatMul(MulTable()), ScalarMatMul(LogTable16()), TableMatMul(),
                           CLMulMatMul())

let rng = Xoshiro(2027)
    for (name, F) in BULK_LEVELS
        g = SUITE["bulk " * name] = BenchmarkGroup()
        GROUP_SIZE["bulk " * name] = N_BULK
        a, b = rand(rng, F, N_BULK), rand(rng, F, N_BULK)
        c, k = similar(a), rand(rng, F)
        m = MulBy(k)
        g["rand!"] = @benchmarkable rand!($rng, $c)
        g["constant * array", "c .* a"] = @benchmarkable broadcast!(*, $c, $k, $a)
        g["constant * array", "MulBy(c).(a)"] = @benchmarkable broadcast!($m, $c, $a)
        g["constant * array", "MulBy(c) setup"] = @benchmarkable MulBy($k)
        if F in BITSLICED_LEVELS
            A, B = BitslicedVector(a), BitslicedVector(b)
            C = similar(A)
            g["bitsliced", "mul_elementwise!"] = @benchmarkable mul_elementwise!($C, $A, $B)
            g["bitsliced", "pack"] = @benchmarkable BitslicedVector($a)
            g["bitsliced", "unpack"] = @benchmarkable Vector($A)
            g["bitsliced", "add"] = @benchmarkable $A .+ $B
        end

        g = SUITE["matmul " * name] = BenchmarkGroup()
        GROUP_SIZE["matmul " * name] = MATMUL_N^3
        A, B = rand(rng, F, MATMUL_N, MATMUL_N), rand(rng, F, MATMUL_N, MATMUL_N)
        C = similar(A)
        for s in MATMUL_STRATEGIES
            sname = s isa ScalarMatMul ? "ScalarMatMul(" * string(nameof(typeof(s.strategy))) * ")" :
                    string(nameof(typeof(s)))
            g[sname] = @benchmarkable matmul!($s, $C, $A, $B)
        end
        g["default (*)"] = @benchmarkable $A * $B
    end
end
