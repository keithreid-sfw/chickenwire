# chickenwire
# copyright Dr Keith Reid Cailleach Computing Ltd 2022

# in tty get REPL
# julia> include("chickenwire.jl")
# julia> p=chickenwire()
# gui()


#s0 libraries
using Plots
using Random
using Test

#s1 config
function get_chw_size()
    chw_size::Int64    = 64 
    return chw_size
end 

function get_chw_n()
    chw_n::Int64       = 8
    return chw_n
end

#s2 model

function build_points(chw_n::Int64)
    chw_size::Int64                 	= get_chw_size() 
    xs::Array{Float64,1}            	= rand(Float64,chw_n)
    ys::Array{Float64,1}            	= rand(Float64,chw_n)
    big_xs::Array{Float64,1}            = [x*chw_size for x in xs] 
    big_ys::Array{Float64,1}            = [y*chw_size for y in ys] 
    points::Array{Array{Float64,1},1}	= [[x[1],x[2]] for x in zip(big_xs, big_ys)] 
    return points
end

#s3 view

function draw_chickenwire(   oranges::Array{Array{Float64,1},1},	 
		       blues::Array{Array{Float64,1},1}, 
		       index::Array{Array{Float64,1},1})
    p=scatter([xy[1] for xy in oranges],[xy[2] for xy in oranges],color="orange",legend=false)
    scatter!([xy[1] for xy in blues],    [xy[2] for xy in blues], color="blue")
    scatter!([xy[1] for xy in index],    [xy[2] for xy in index], color="red")
    return p
end

function find_tether(this_point::Array{Float64,1})
    point_x                 = this_point[1]
    point_y                 = this_point[2]
    int_x::Int64            = trunc(Int64, point_x)
    int_y::Int64            = trunc(Int64, point_y)
    residual_x              = point_x - int_x
    residual_y              = point_y - int_y
    int_int::Array{Int64,1} = [int_x, int_y]
    tether                  = int_int
    subcell_type = "top"
    if subcell_type == "top" 
        tether = int_int
    end
    return tether
end

function find_tethers(some_points::Array{Array{Float64,1},1})
    tethers::Array{Array{Int64,1},1} = []
    for this_point in some_points
        this_tether::Array{Int64,1} = find_tether(this_point)
        push!(tethers, this_tether)
    end
    return tethers 
end

function find_bin_ends(point::Array{Float64,1})
    point_x     		= trunc(Int64, point[1])
    x_as_2_bin			= bitstring(point_x)[length(bitstring(point_x))-1:end]
    point_y     		= trunc(Int64, point[2])
    y_as_2_bin	      		= bitstring(point_y)[length(bitstring(point_y))-1:end]
    bin_ends::Array{String,1}   = [x_as_2_bin, y_as_2_bin]
    return bin_ends
end

function find_subcell_shape(bin_ends::Array{String,1})
    x_end = bin_ends[1]
    println("first bitstring:\t", x_end)
    if x_end[2] == 1
        print("odd!")
        shape = "s" # its a square
    else
        print("even!")
        shape = "t" # its trianglier
    end
    return shape
end

function test_find_subcell_shape()
    
    #=
    8
    7/T\B/T\B
    6\B/T\B/T
    5/T\B/T\B 
    4/T\B/T\B
    3\B/T\B/T
    2/T\B/T\B
    1\B/T\B/T
    0/T\B/T\B
     012345678
    =#

    notional_top        = [1.5, 0.5]
    bin_ends 		= find_bin_ends(notional_top)
    subcell_shape       = find_subcell_shape(bin_ends)
    @test subcell_shape == "s"

    notional_NER        = [0.8, 0.2]
    bin_ends 		= find_bin_ends(notional_NER)
    subcell_shape        = find_subcell_shape(bin_ends)
    @test subcell_shape == "t"

    notional_NEL        = [0.2, 0.8]
    bin_ends 		= find_bin_ends(notional_NEL)
    subcell_shape        = find_subcell_shape(bin_ends)
    @test subcell_shape == "t"

    notional_bottom     = [1.5, 1.5]
    bin_ends 		= find_bin_ends(notional_bottom)
    subcell_shape        = find_subcell_shape(bin_ends)
    @test subcell_shape == "s"

    notional_NWL        = [0.2, 1.2]
    bin_ends 		= find_bin_ends(notional_NWL)
    subcell_shape        = find_subcell_shape(bin_ends)
    @test subcell_shape == "t"

    notional_NWR        = [0.8, 1.8] 
    bin_ends 		= find_bin_ends(notional_NWR)
    subcell_shape        = find_subcell_shape(bin_ends)
    @test subcell_shape == "t"


    println("bin ends okay for six prototypes")

end

    
function test_find_tethers()
    reference = [[1.1,0.2],[1.1,0.2],[1.1,0.2],[1.1,0.2],[1.1,0.2]]
    tethers   = find_tethers(reference)
    @test typeof(tethers) == Array{Array{Int64,1},1} 
    println("tethers type is right")
end

function test_find_bin_ends()
    
    #=
    8
    7/T\B/T\B
    6\B/T\B/T
    5/T\B/T\B 
    4/T\B/T\B
    3\B/T\B/T
    2/T\B/T\B
    1\B/T\B/T
    0/T\B/T\B
     012345678
    =#

    notional_top        = [1.5, 0.5]
    bin_ends 		= find_bin_ends(notional_top)
    @test bin_ends      == ["01","00"]

    notional_NER        = [0.8, 0.2]
    bin_ends 		= find_bin_ends(notional_NER)
    @test bin_ends      == ["00","00"]

    notional_NEL        = [0.2, 0.8]
    bin_ends 		= find_bin_ends(notional_NEL)
    @test bin_ends      == ["00","00"]

    notional_bottom     = [1.5, 1.5]
    bin_ends 		= find_bin_ends(notional_bottom)
    @test bin_ends      == ["01","01"]

    notional_NWL        = [0.2, 1.2]
    bin_ends 		= find_bin_ends(notional_NWL)
    @test bin_ends      == ["00","01"]

    notional_NWR        = [0.8, 1.8] 
    bin_ends 		= find_bin_ends(notional_NWR)
    @test bin_ends      == ["00","01"]
    println("bin ends okay for six prototypes")
end

function all_chickenwire_tests()
    test_find_tethers()
    test_find_bin_ends()
    test_find_subcell_shape()
end

#s4 controller

function chickenwire()
    chw_n::Int64		        = get_chw_n()
    oranges::Array{Array{Float64,1},1} 	= build_points(chw_n)
    blues::Array{Array{Float64,1},1}    = build_points(chw_n)
    index::Array{Array{Float64,1},1}	= build_points(1::Int64)
    
    orange_bin_ends = [find_bin_ends(orange) for orange::Array{Float64} in oranges] 
    blue_bin_ends = [find_bin_ends(blue) for blue::Array{Float64} in blues] 
    index_bin_ends = [find_bin_ends(singleton) for singleton::Array{Float64} in index] 





    p = draw_chickenwire(oranges,
			 blues,   
		         index)
    orange_tethers = find_tethers(oranges)
    p = scatter!([xy[1] for xy in orange_tethers], [xy[2] for xy in orange_tethers])
    println("oranges:\n", oranges)
    println("blues:\n",   blues)
    println("index:\n",   index)
    println("orange bin ends:\n", orange_bin_ends)
    println("blues bin ends:\n",   blue_bin_ends)
    println("index bin ends:\n",   index_bin_ends)
    all_chickenwire_tests()
    return(p)
end
