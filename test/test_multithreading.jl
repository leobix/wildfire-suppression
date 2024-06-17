using Profile

# test multithreading
Threads.@threads for i ∈ 1:100
	println(Threads.threadid())
	j = 0.0
	for k ∈ 1:10000
		for l ∈ 1:1000
			j += (l + k)^0.5
		end
	end
end

println()
println("again")
println()

Profile.init()
@profile begin
	Threads.@threads for i ∈ 1:100
		println(Threads.threadid())
		j = 0.0
		for k ∈ 1:10000
			for l ∈ 1:1000
				j += (l + k)^0.5
			end
		end
	end
end
Profile.print(mincount = 5)
