macro timeit(ex)
    return quote
        local t0 = time()
        local val = $(esc(ex))
        println("Elapsed time: ", time() - t0, " seconds")
        val
    end
end