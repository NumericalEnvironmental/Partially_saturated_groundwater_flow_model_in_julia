### input and output functions for jFlow.jl


using DelimitedFiles


function ReadMaterials()::Array{Material, 1}
    # read various model parameters from file
    material = Material[]
    data = readdlm("Materials.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        index = Int64(data[1][i, 1])
        Kh = Float64(data[1][i, 2])                 # horizontal and vertical hydraulic conductivities
        Kz = Float64(data[1][i, 3])
        Ss = Float64(data[1][i, 4])                 # specific storage
        phi = Float64(data[1][i, 5])                # porosity
        alpha = Float64(data[1][i, 6])              # Van Genuchten unsaturated hydraulic parameters
        N = Float64(data[1][i, 7])
        m = 1.0 - 1.0/N
        push!(material, Material(index, Kh, Kz, Ss, phi, alpha, N, m))
    end
    println("Read material properties.")
    return material
end


function ReadParams()::Params
    # read numerical model parameters
    data = readdlm("Params.txt", '\t', header=false)
    endTime = Float64(data[1, 2])                   # end of simulation
    dt0 = Float64(data[2, 2])                       # starting time step
    dtMax = Float64(data[3, 2])                     # maximum time step
    dhMax = Float64(data[4, 2])                     # maximum head change (anywhere) per time step
    f = Float64(data[5, 2])                         # time step adjustment factor
    c = Float64(data[6, 2])                         # explicit time-step weighting (0=explicit, 0.5=Crank-Nicholson, 1=fully-implicit)
    params = Params(endTime, dt0, dtMax, dhMax, f, c)
    println("Read model run parameters.")
    return params
end


function ReadCells(material::Array{Material, 1})::Array{Cell, 1}
    # read cell setup file and process cell objects
    cell = Cell[]
    data = readdlm("Cells.txt", '\t', header=true)
    for iRow = 1:size(data[1], 1)
        cellIndex = Int64(data[1][iRow, 1])
        x0 = Float64(data[1][iRow, 2])              # locaiton of volume lement/cell node point
        y0 = Float64(data[1][iRow, 3])
        z0 = Float64(data[1][iRow, 4])
        vol = Float64(data[1][iRow, 5])             # cell volume; note: set equal to very large number for fixed-head condition
        b = Float64(data[1][iRow, 6])               # cell vertical thickness (currently not in use)
        h = Float64(data[1][iRow, 7])               # hydraulic head (an initial condition, at start up)
        matIndex = Int64(data[1][iRow, 8])          # index number for material
        mat = material[matIndex]
        xNum = Int64(data[1][iRow, 9])              # number of identical cells along different axes    
        yNum = Int64(data[1][iRow, 10])     
        zNum = Int64(data[1][iRow, 11])
        xInc = Int64(data[1][iRow, 12])             # increment number for cellIndex while populating along different axes  
        yInc = Int64(data[1][iRow, 13])     
        zInc = Int64(data[1][iRow, 14])     
        xStep = Float64(data[1][iRow, 15])          # change in node point coordinate while populating along different axes 
        yStep = Float64(data[1][iRow, 16])      
        zStep = Float64(data[1][iRow, 17])
        for k = 1:zNum, j = 1:yNum, i = 1:xNum
            index = cellIndex + (i-1)*xInc + (j-1)*yInc + (k-1)*zInc
            x = x0 + (i-1)*xStep
            y = y0 + (j-1)*yStep
            z = z0 + (k-1)*zStep
            connects = Int64[]
            neighbors = Int64[]
            push!(cell, Cell(index, x, y, z, vol, b, h-z, 0.0, mat, connects, neighbors))
            cell[end].S = Sat(cell[end])        # calculate initial saturation, now that cell has been defined
        end
    end
    println("Read and assigned cell properties.")
    return cell
end


function ReadCellMods(cell::Array{Cell, 1}, material::Array{Material, 1})::Array{Cell, 1}
    # modify select cell properties, if specified
    data = readdlm("CellMods.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        cellIndex = Int64(data[1][i, 1]) 
        vol = Float64(data[1][i, 2])
        h0 = Float64(data[1][i, 3])
        matIndex = Int64(data[1][i, 4])
        cell[cellIndex].vol = vol
        cell[cellIndex].P = h0 - cell[cellIndex].z
        cell[cellIndex].mat = material[matIndex]
    end
    println("Modified subset of cells, as specified.")  
    return cell
end


function ReadGridConnects()::Array{Connection, 1}
    # read and assign connections, formatted for internal grid cells
    connect = Connection[]
    data = readdlm("gridConnections.txt", '\t', header=true)
    for iRow = 1:size(data[1], 1)
        cellStart = Int64(data[1][iRow, 1])                 # starting cell index for connection sequence
        numCellStart = Int64(data[1][iRow, 2])              # number of start cells (along row, col, etc.)
        stepCellStart = Int64(data[1][iRow, 3])             # increment between start cells     
        numPopulate = Int64(data[1][iRow, 4])               # number of connections in direction of population
        stepPopulate = Int64(data[1][iRow, 5])              # increment between populated cells 
        del1 = Float64(data[1][iRow, 6])                    # node point-to-interface distances and interface area
        del2 = Float64(data[1][iRow, 7])
        A = Float64(data[1][iRow, 8])       
        zConnect = Bool(data[1][iRow, 9])                   # boolean flag indicating if connect is oriented vertically
        for j = 1:numCellStart, i = 1:numPopulate
            cell1 = cellStart + (j-1)*stepCellStart + (i-1)*stepPopulate
            cell2 = cellStart + (j-1)*stepCellStart + i*stepPopulate
            push!(connect, Connection(cell1, cell2, del1, del2, A, zConnect, 0.0))           
        end
    end
    println("Read and assigned internal grid cell-to-cell connections.")
    return connect
end


function ReadSpecConnects(connect::Array{Connection, 1})::Array{Connection, 1}
    # read and assign special cell connections for boundaries, unique cells, etc.
    data = readdlm("specConnections.txt", '\t', header=true)
    for iRow = 1:size(data[1], 1)
        cellStart = Int64(data[1][iRow, 1])                 # starting cell index for connection sequence
        numCellStart = Int64(data[1][iRow, 2])              # number of start cells (along row, col, etc.)
        stepCellStart = Int64(data[1][iRow, 3])             # increment between start cells     
        boundaryCell = Int64(data[1][iRow, 4])              # index number of cell completing the connection
        del1 = Float64(data[1][iRow, 5])
        del2 = Float64(data[1][iRow, 6])
        A = Float64(data[1][iRow, 7])       
        zConnect = Bool(data[1][iRow, 8])       
        for j = 1:numCellStart
            cell1 = cellStart + (j-1)*stepCellStart
            cell2 = boundaryCell
            push!(connect, Connection(cell1, cell2, del1, del2, A, zConnect, 0.0))           
        end
    end
    println("Read and assigned special cell-to-cell connections.")
    return connect
end


function ReadMonitors()::Array{Int64, 1}
    # return a list of monitoring point (cell) indices
    data = readdlm("Monitors.txt", '\t', header=true)
    monitor = Int.(data[1])
    monitor = vec(monitor)
    return monitor
end


function ReadSources()
    # read source/sinks terms and assign to cells & stress periods
    data = readdlm("Sources.txt", '\t', header=true)
    stressPeriods = sort(collect(Set(float(data[1][:, 2]))))
    sourceCells = sort(collect(Set(data[1][:, 1])))
    sourceCells = floor.(Int, sourceCells)
    source = zeros(Float64, length(sourceCells), length(stressPeriods))
    for i = 1:size(data[1], 1)
        cellIndex = Int64(data[1][i, 1]) 
        t0 = Float64(data[1][i, 2])             # starting time for source
        Q = Float64(data[1][i, 3])              # volumetric flux (L^3/T)
        iCol = FindStressPeriod(t0, stressPeriods)
        iRow = findall(sourceCells.==cellIndex)[1]
        source[iRow, iCol:end].= Q
    end
    println("Processed source terms.")
    return stressPeriods, sourceCells, source
end


function WriteCells(cell::Array{Cell, 1}, fileName::String, boundFlag::Int64)
    # summarize cell properties to file
    csvfile = open(fileName,"w")
    line_out = "cellIndex" * "," * "x" * "," * "y" * "," * "z" * "," *
        "vol" * "," * "P" * "," * "S" * "," * "material"
    println(csvfile, line_out)
    for (i, ce) in enumerate(cell)
        if (ce.vol < 1e+10) || boundFlag == 1       # if applicable, exclude boundary elements from summary
            ce.S = Sat(ce)      # update saturation estimate
            line_out = string(ce.index) * "," * string(ce.x) * "," * string(ce.y) * "," * string(ce.z) * "," *
                string(ce.vol) * "," * string(ce.P) * "," * string(ce.S) * "," * string(ce.mat.index)
            println(csvfile, line_out)
        end
    end
    close(csvfile)
    println("Wrote " * fileName * ".")
end


function WriteConnections(connect::Array{Connection, 1})
    # summarize cell properties to file
    csvfile = open("ConnectionSummary.csv","w")
    line_out = "connection" * "," * "cell1" * "," * "cell2" * "," * "del1" * "," *
        "del2" * "," * "area" * "," * "zConnect"
    println(csvfile, line_out)
    for (i, cn) in enumerate(connect)
        line_out = string(i) * "," * string(cn.cell1) * "," * string(cn.cell2) * "," * string(cn.del1) * "," * 
            string(cn.del2) * "," * string(cn.A) * "," * string(cn.zConnect)
        println(csvfile, line_out)
    end
    close(csvfile)
    println("Wrote connections file.")
end


function WriteSources(stressPeriod::Array{Float64, 1}, sourceCells::Array{Int64, 1}, source::Array{Float64, 2})
    # summarize source term(s) to file
    csvfile = open("SourceSummary.csv","w")
    line_out = "cell"
    for t in stressPeriod
        line_out *= "," * string(t)
    end
    println(csvfile, line_out)  
    for cellIndex in sourceCells
        line_out = string(cellIndex)
        iRow = findall(sourceCells.==cellIndex)[1]
        for iCol = 1:length(stressPeriod)
            line_out *= "," * string(source[iRow, iCol])
        end
        println(csvfile, line_out)
    end
    close(csvfile)
    println("Wrote source terms file.")
end


function WriteTimeSeries(monitor::Array{Int64, 1}, timeMon::Array{Float64, 1}, timeSeries)
    # write monitoring point time series to file
    csvfile = open("TimeSeries.csv","w")
    line_out = "time"
    for m in monitor
        line_out *= "," * "cell" * "_" * string(m)
    end
    println(csvfile, line_out)  
    for (i, t) in enumerate(timeMon)
        line_out = string(t)
        for j = 1:length(monitor)
            line_out *= "," * string(timeSeries[j][i])
        end
        println(csvfile, line_out)
    end
    close(csvfile)
    println("Wrote time series file.")    
end


function WriteVectors(cell::Array{Cell, 1}, vx::Array{Float64, 1}, vy::Array{Float64, 1}, vz::Array{Float64, 1})
    # write velocity vectors at cell node points
    csvfile = open("Velocity.csv","w")
    line_out = "x" * "," * "y" * "," * "z" * "," *
        "vx" * "," * "vy" * "," * "vz"
    println(csvfile, line_out)
    for (i, ce) in enumerate(cell)
        if (ce.vol < 1e+10)         # exclude boundary elements
            line_out = string(ce.x) * "," * string(ce.y) * "," * string(ce.z) * "," *
                string(vx[i]) * "," * string(vy[i]) * "," * string(vz[i])
            println(csvfile, line_out)
        end
    end
    close(csvfile)
    println("Wrote velocity vectors file.")
end

