const a = [984,981,976,950,899,890,887,885,880,800,798,790,777,767,750,701,
    697,688,680,678,650,599,589,567,550,501,9,8,7,6,5,4,3,2,1]
const b = [1090,1080,1074,1065,1057,1056,1047,1041,1041,1038,1025,1013,1008,
    992,991,991,991,978,977,959,945,935,925,925,923,915,908,904,901,901,
    900,897,894,882,880,876,866,854,849,849,833,818,818,812,811,809,798,
    794,793,788,772,763,747,746,743,737,736,734,732,730,728,718,714,713,
    706,701,699,691,690,689,681,672,663,656,654,653,652,651,646,644,640,
    637,637,635,634,633,630,625,621]

struct DigitSet
    x::UInt32
end

function DigitSet(x::Integer)
    set = DigitSet()
    absolute = unsigned(abs(x))
    while absolute > 0
        set = push!(set, rem(absolute, 10))
        absolute = oftype(absolute, div(absolute, 10))
    end
    return set
end

DigitSet() = DigitSet(zero(UInt32))
Base.push!(x::DigitSet, y::Integer) = DigitSet(x.x | (UInt32(1) << y))
Base.union(x::DigitSet, y::DigitSet) = DigitSet(x.x | y.x)
Base.hash(x::DigitSet, h::UInt) = hash(x.x, h)
isdisjoint(x::DigitSet, y::DigitSet) = x.x & y.x == zero(UInt32)
Base.issubset(x::DigitSet, y::DigitSet) = (x.x & ~y.x) == zero(UInt32)

function linear_deduplicate(numbers)
    bydigits = Dict{DigitSet, eltype(numbers)}()
    for number in numbers
        set = DigitSet(number)
        existing = get(bydigits, set, typemin(Int))
        if number > existing
            bydigits[set] = number
        end
    end
    return collect(bydigits)
end

function square_deduplicate(numberpairs)
    accepted = Int[]
    for (set, number) in numberpairs
        add = true
        for (set_, number_) in numberpairs
            if issubset(set_, set) && number_ > number
                add = false
                break
            end
        end
        if add
            push!(accepted, number)
        end
    end
    # We sort to prune the tree of recusive possibilities earlier
    sort!(accepted, rev=true)
    return accepted
end

abstract type Solution end

struct SmallSolution <: Solution
    sum::Int64
    x::UInt64
end

mutable struct BigSolution <: Solution
    sum::Int64
    x::Vector{Int}
end

SmallSolution(x) = SmallSolution(x, x)
SmallSolution() = SmallSolution(0)
BigSolution() = BigSolution(0, Int[])

Base.copy(x::SmallSolution) = x
Base.copy(x::BigSolution) = BigSolution(x.sum, copy(x.x))

function pushindex!(x::SmallSolution, numbers, index)
    @inbounds n = numbers[index]
    return SmallSolution(n+x.sum, x.x | UInt64(1) << (index - 1))
end

function pushindex!(x::BigSolution, numbers, index)
    n = numbers[index]
    x.sum += n
    push!(x.x, n)
    return x
end

Base.collect(x::BigSolution, nums::Vector{T}) where T<:Integer = x.x

function Base.collect(x::SmallSolution, nums::Vector{T}) where T<:Integer
    result = Int[]
    @inbounds for i in 1:min(64, unsigned(length(nums)))
        if x.x >>> (i-1) & UInt64(1) == UInt64(1)
            push!(result, nums[i])
        end
    end
    return result
end

function _recurse(current::Solution, numbers, sets, index, set, solution)
    if index > length(numbers)
        if current.sum > solution[].sum
            solution[] = current
        end
        return nothing
    end
    @inbounds nextset = sets[index]
    _recurse(current, numbers, sets, index+1, set, solution)
    if isdisjoint(set, nextset)
        nextcurrent = pushindex!(copy(current), numbers, index)
        _recurse(nextcurrent, numbers, sets, index+1, union(set, nextset), solution)
    end
    return nothing
end

function recurse(numbers)
    numbersarr = collect(numbers)
    T = length(numbersarr) ≤ 64 ? SmallSolution : BigSolution
    solution = Ref(T())
    _recurse(T(), numbersarr, map(DigitSet, numbersarr), 1, DigitSet(), solution)
    return sort!(collect(solution[], numbers))
end

function solve(numbers::Vector{Int})
    filtered = linear_deduplicate(numbers)
    filtered2 = square_deduplicate(filtered)
    return recurse(filtered2)
end
