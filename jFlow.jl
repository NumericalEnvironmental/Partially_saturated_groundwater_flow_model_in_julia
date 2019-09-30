########################################################################################
#
# jFlow: a 3-D saturated-unsaturated groundwater flow model
#
# (1) Methodology is the integral finite difference method (IFDM)
# (2) Capillary forces are ignored (i.e., intended for larger spatial scales)
# (3) Pressure head in a cell is defined at z of node point (i.e., if P = 0, S is typically = 0.5)
#
########################################################################################


using SparseArrays


### data types


mutable struct Material
    index::Int64
    Kh::Float64                         # saturated hydraulic conductivity components
    Kz::Float64
    Ss::Float64                         # specific storage
    phi::Float64                        # porosity
    alpha::Float64                      # Van Genuchten model parameters
    N::Float64
    m::Float64
end


mutable struct Cell
    index::Int64                        # index tracking number for cell (for initial struct array setup)
    x::Float64                          # grid cell center coordinates
    y::Float64
    z::Float64                          # node point elevation
    vol::Float64                        # cell volume
    b::Float64                          # cell vertical thickness
    P::Float64                          # pressure head
    S::Float64                          # saturation
    mat::Material                       # corresponding material
    connects::Array{Int64, 1}           # list of associated connection indices
    neighbors::Array{Int64, 1}          # list of connected cells
end


mutable struct Connection
    cell1::Int64                        # connected cells (index numbers)
    cell2::Int64
    del1::Float64                       # respective distances from cell centerpoints to interface
    del2::Float64
    A::Float64                          # interfacial area of the connection
    zConnect::Bool                      # connection is designated as vertically-oriented
	KMean::Float64 						# mean saturated hydraulic conductivity across connection
end


mutable struct Params
    endTime::Float64                    # end of simulation time
    dt0::Float64                        # initial (and minimum) time step
    dtMax::Float64                      # maximum time step
    dhMax::Float64                      # maximum head change in any model cell, per time step
    f::Float64                          # time step multiplier
    c::Float64                          # Crank-Nicolson implicit solution weighting factor
end


include("jFlowIOModule.jl")             # various input and output functions included in separate module


### hydraulic parameterization functions


function kr(cell::Cell)::Float64
    # relative permeability as a function of water saturation
    return cell.S^0.5 * (1.0 - (1.0 - cell.S^(1.0/cell.mat.m))^cell.mat.m)^2
end


function Sat(cell::Cell)::Float64
    # (effective) saturation as a function of pressure head
    if cell.P < 0.0
        S = (1. + abs(cell.mat.alpha * cell.P)^cell.mat.N)^(-cell.mat.m)          
    else
        S = 1.0
    end
    return S
end


function dSdP(cell::Cell)::Float64
    # change in effective saturation per change in capillary pressure head
    dS = -(cell.mat.alpha*abs(cell.P))^cell.mat.N * ((cell.mat.alpha*abs(cell.P))^cell.mat.N + 1.0)^(-cell.mat.m-1.0) * cell.mat.m * cell.mat.N/cell.P
    return dS * cell.mat.phi
end


function Storage(cell::Cell)::Float64
    # return appropriate storage term, depending on saturation state of cell
    if Sat(cell) < 1.0
        stor = dSdP(cell) * cell.vol
    else
        stor = cell.mat.Ss * cell.vol
    end
    return stor
end


### utility functions


function SortCell(tempCell::Array{Cell, 1})::Array{Cell, 1}
    # insertion sort approach to create ordered array of cell structs
    cell = [tempCell[1]]
    seq = [tempCell[1].index]
    for i = 2:length(tempCell)
        index = length(findall(seq.<tempCell[i].index)) + 1
        insert!(seq, index, tempCell[i].index)
        insert!(cell, index, tempCell[i])
    end
    println("Sorted cell struct array.")
    return cell
end


function MatchConnects(cell::Array{Cell, 1}, connect::Array{Connection, 1})::Array{Cell, 1}
    # populate cell connection lists
    for (i, cn) in enumerate(connect)
        push!(cell[cn.cell1].connects, i)
        push!(cell[cn.cell2].connects, i)
        push!(cell[cn.cell1].neighbors, cn.cell2)
        push!(cell[cn.cell2].neighbors, cn.cell1)       
    end
    return cell
end


function FindStressPeriod(t::Float64, stressPeriods::Array{Float64, 1})::Int64
    # index number of stress period to which "t" belongs
    return sum((t.>=stressPeriods)*1)
end


function BlendedK(connect::Array{Connection, 1}, cell::Array{Cell, 1})::Array{Connection, 1}
    # weighted harmonic mean for saturated hydraulic conductivity heterogeneous cell interfaces
	for cn in connect
		if cn.zConnect == false
			cn.KMean = (cn.del1 + cn.del2)/(cn.del1/cell[cn.cell1].mat.Kh + cn.del2/cell[cn.cell2].mat.Kh)
		else
			cn.KMean = (cn.del1 + cn.del2)/(cn.del1/cell[cn.cell1].mat.Kz + cn.del2/cell[cn.cell2].mat.Kz)    
		end
	end
    return connect
end


function EffectiveK(cn::Connection, cell::Array{Cell, 1})::Float64
    # combo weighted harmonic mean & arithmetic means for unsaturated cell interfaces
    krMean = (cn.del1*kr(cell[cn.cell1]) + cn.del2*kr(cell[cn.cell2]))/(cn.del1 + cn.del2)
    return krMean * cn.KMean
end


function TimeStep(t::Float64, dt::Float64, params::Params, tStops::Array{Float64, 1})
    # revise time step
    finFlag = false
    dtNew = minimum([dt*params.f, params.dtMax])    
    for tStop in tStops
        if (t < tStop) && (t+dtNew > tStop)
            finFlag = (t+dtNew > tStops[end])       # indicates end-of-simulation-time has been reached      
            dtNew = tStop - t
            break
        end
    end
    return dtNew, finFlag
end


function Vel(cell::Array{Cell, 1}, connect::Array{Connection, 1})
    # compute pore velocity vector arrays at cell node points
    vx = zeros(Float64, length(cell))
    vy = zeros(Float64, length(cell))    
    vz = zeros(Float64, length(cell))
    xCompCell = zeros(Float64, length(cell))    # connection component projections along respective axes, per cell
    yCompCell = zeros(Float64, length(cell))
    zCompCell = zeros(Float64, length(cell))    
    for cn in connect
        conduct = EffectiveK(cn, cell)        # effective hydraulic conductivity
        d = cn.del1+cn.del2
        xComp = (cell[cn.cell2].x - cell[cn.cell1].x)/d^2
        yComp = (cell[cn.cell2].y - cell[cn.cell1].y)/d^2        
        zComp = (cell[cn.cell2].z - cell[cn.cell1].z)/d^2        
        q = conduct * ((cell[cn.cell1].P+cell[cn.cell1].z)-(cell[cn.cell2].P+cell[cn.cell2].z)) / d
        vx[cn.cell1] += xComp * q/(cell[cn.cell1].mat.phi*cell[cn.cell1].S) 
        vx[cn.cell2] += xComp * q/(cell[cn.cell2].mat.phi*cell[cn.cell2].S)        
        vy[cn.cell1] += yComp * q/(cell[cn.cell1].mat.phi*cell[cn.cell1].S)
        vy[cn.cell2] += yComp * q/(cell[cn.cell2].mat.phi*cell[cn.cell2].S)     
        vz[cn.cell1] += zComp * q/(cell[cn.cell1].mat.phi*cell[cn.cell1].S)
        vz[cn.cell2] += zComp * q/(cell[cn.cell2].mat.phi*cell[cn.cell2].S)
        xCompCell[cn.cell1] += abs(xComp)
        xCompCell[cn.cell2] += abs(xComp)
        yCompCell[cn.cell1] += abs(yComp)
        yCompCell[cn.cell2] += abs(yComp)
        zCompCell[cn.cell1] += abs(zComp)
        zCompCell[cn.cell2] += abs(zComp)        
    end
    # normalize velocity components at each cell node by number of connections
    for i = 1:length(cell)
        vx[i] /= (xCompCell[i] + 1e-10)                 # 1e-10 factor added to prevent divide-by-zero errors
        vy[i] /= (yCompCell[i] + 1e-10)
        vz[i] /= (zCompCell[i] + 1e-10)        
    end
    println("Computed pore velocity vectors.")
    return vx, vy, vz
end


### matrix operations functions


function AssembleMatrix(cell::Array{Cell, 1}, connect::Array{Connection, 1}, params::Params, dt::Float64, t::Float64, stressPeriods::Array{Float64, 1}, sourceCells::Array{Int64, 1}, source::Array{Float64, 2})
    # fill out the LHS of the equation matrix
    row_index = Int64[]                     # indexing system for sparse matrix
    col_index = Int64[]
    data = Float64[]
    b = Float64[]
    for (i, ce) in enumerate(cell)          # for each row
        diag = Storage(ce)/dt
        bSum = 0.
        for (j, nghbr) in enumerate(ce.neighbors)           # for each neighboring cell
            conduct = EffectiveK(connect[ce.connects[j]], cell) *
                connect[ce.connects[j]].A / (connect[ce.connects[j]].del1 + connect[ce.connects[j]].del2)
            push!(row_index, i)                                             # left-hand-side matrix
            push!(col_index, nghbr)
            push!(data, -params.c * conduct)
            diag += params.c * conduct
            bSum += conduct * ((cell[nghbr].P + cell[nghbr].z) - (ce.P + ce.z))         # right-hand-side (explicit) vector
        end
        push!(row_index, i)         # set diagonal term
        push!(col_index, i)
        push!(data, diag)        
        push!(b, bSum)              # set explicit vector
    end
    # add source/sink terms to b[], as warranted
    j = FindStressPeriod(t, stressPeriods)
    for (i, cellIndex) in enumerate(sourceCells)
        b[cellIndex] += source[i, j]
    end
    return data, row_index, col_index, b
end
    

### main script


function jFlow(mode::Int64)

    # read and process input files
    material = ReadMaterials()                          # read materials properties
    tempCell = ReadCells(material)                      # read and assign cells (possibly out of sequence)
    cell = SortCell(tempCell)                           # use insertion sort approach to sequence cell struct array
    #cell = ReadCellMods(cell, material)                # modify specific cell properties, as specified
    connect = ReadGridConnects()                        # read internal grid connections
    #connect = ReadSpecConnects(connect)                 # append unique or boundary connections
    cell = MatchConnects(cell, connect)                 # assign cell connection lists
	connect = BlendedK(connect, cell) 					# pre-calculated connection saturated K's
    stressPeriods, sourceCells, source = ReadSources()  # tabulate stress periods and source/sink cells
	drainCell, drainP = ReadDrains() 					# read cells designated as drains
    params = ReadParams()                               # read model properties
    monitor = ReadMonitors()                            # note monitor locations (results tabulated for every time step)
    
    # create container arrays for monitoring points
    timeMon = Float64[]
    timeSeries = []
    for m in monitor
        push!(timeSeries, Float64[]) 
    end
    
    # summarize model setup
    WriteCells(cell, "CellSummary.csv", 1, drainCell, drainP)
    WriteConnections(connect)
    WriteSources(stressPeriods, sourceCells, source)
   

   
    if mode != 0            # if mode == 0, do not run model (just print out model setup and exit)
    
        # initialize model
        t = 0.0
        dt = params.dt0                     # first time step size
        tStops = copy(stressPeriods)        # times for which dt must be adjusted for "t" to match ...
        push!(tStops, params.endTime)
        finFlag = false
        convgFlag = false
        dP = zeros(Float64, length(cell))
        
        while finFlag==false
        
            while convgFlag == false
            
                # assemble matrix
                data, row_index, col_index, b = AssembleMatrix(cell, connect, params, dt, t, stressPeriods, sourceCells, source)
                A = sparse(row_index, col_index, data, length(cell), length(cell))
                
                # solve system of equations
                dP = \(A, b)
                
                # adjust cells with drains
				for (i, indx) in enumerate(drainCell)
					dP[indx] = minimum([dP[indx], drainP[i]-cell[indx].P])
				end
                
                # check convergence (maximum head change); adjust dt if needed
                if maximum(abs.(dP)) >= params.dhMax
                    dt /= params.f
                    @assert(dt > params.dt0)             # stop execution if dt falls below minimum
                else
                    convgFlag = true
                end
                
            end
            
            # update results
            t += dt
            for (i, ce) in enumerate(cell)
                ce.P += dP[i]
                ce.S = Sat(ce)
            end
            println("\t", "t = ", string(t))
            
            # write to monitor file(s)
            push!(timeMon, t)
            for (i, monCell) in enumerate(monitor)
                push!(timeSeries[i], cell[monCell].P)
            end
            
            # determine next time step
            dt, finFlag = TimeStep(t, dt, params, tStops)
        
            convgFlag = false           # reset head-change flag for next time step
        
        end
    
    end

    # write output
    WriteCells(cell, "FinalState.csv", 0, drainCell, drainP)	# write cell output file
    vx, vy, vz = Vel(cell, connect)                 			# compute and write velocity vectors to file
    WriteVectors(cell, vx, vy, vz)
    WriteTimeSeries(monitor, timeMon, timeSeries)   			# write monitor point time series to file
    
    println("Finished.")
    
end

### run script
jFlow(1)