# test multithreading
Threads.@threads for i ∈ 1:100
    println(Threads.threadid())
    j = 0
    for k ∈ 1:10000000
        j += 1
    end
end

println()
println("again")
println()

Threads.@threads for i ∈ 1:100
    println(Threads.threadid())
    j = 0
    for k ∈ 1:10000000
        j += 1
    end
end