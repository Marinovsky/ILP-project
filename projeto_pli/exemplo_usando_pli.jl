using JuMP
using Gurobi

function Exemplo_PLI(T::TriggerArcTSP, MaxTime::IntType)
	 global verbose_mode = true
	 global maximum_number_of_nodes_to_print
	(InArcs,OutArcs) = Build_In_Out_Arcs(T)
	AssoArcs = Build_Associated_Arcs(T)

	if (verbose_mode)
	   #Print_Incident_Arcs(T,InArcs,"InArc")
	   #Print_Incident_Arcs(T,OutArcs,"OutArc")
	   #Print_Associated_Arcs(T, AssoArcs, "AssoArcs")
	end

	if (verbose_mode)
   	   PrintHeader("Informacoes da execucao do solver")
	end
	TSP_model = Model(Gurobi.Optimizer)
	if (!verbose_mode)
	   set_optimizer_attribute(TSP_model, "OutputFlag", 0)

	   set_silent(TSP_model)
	end
	
	set_time_limit_sec(TSP_model,MaxTime)
	
	StartingTime = time()

	# x_a
   	@variable(TSP_model, x[i=1:T.NArcs], Bin)
	for a in 1:T.NArcs
	    set_name(x[a], "x_"*string(T.Arc[a].u)*"_"*string(T.Arc[a].v) )
	end

	# y_r
	@variable(TSP_model, y[i=1:T.NTriggers], Bin)
	for r in 1:T.NTriggers
	    set_name(y[r], "y_"*string(T.Trigger[r].trigger_arc_id)*"_"*string(T.Trigger[r].target_arc_id) )
	end

	# z_a
	@variable(TSP_model, z[i=1:T.NArcs], Bin)
	for a in 1:T.NArcs
		set_name(z[a], "z_"*string(T.Arc[a].u)*"_"*string(T.Arc[a].v) )
	end

	# w_r
	@variable(TSP_model, w[i=1:T.NTriggers], Bin)
	for r in 1:T.NTriggers
		set_name(w[r], "w_"*string(T.Trigger[r].trigger_arc_id)*"_"*string(T.Trigger[r].target_arc_id) )
	end

	# u_j
	# The variables u_j represent the respective order of visit for each node in the tour
	@variable(TSP_model, 0 <= u[i=1:T.NNodes] <= T.NNodes, Int)

	# Constraint 1
	# Total sum of x_a is equal to the number of nodes
	expr = 0 # Start an empty expression
	for a in 1:T.NArcs
		expr += x[a]
	end
	@constraint(TSP_model, expr==T.NNodes)

	# Constraint 2
	M = T.NNodes + 1
	println("==============================================")
	println("M: ", M)
	println("==============================================")
	for a in 1:T.NArcs
		if T.Arc[a].v == 1
			continue
		end
		i=T.Arc[a].u
		j=T.Arc[a].v
		@constraint(TSP_model, u[i] - u[j] + M*x[a]<= -1 + M)
	end

   	# Constraint 3
	# In degree constraint: Total sum of x_a, with arcs entering a node, is equal to 1.
   	for v in 1:T.NNodes
	    expr=0 # Start an empty expression
		for a in InArcs[v]
			expr += x[a]
	    end
      	@constraint(TSP_model, expr==1)
    end
	
	# Constraint 4
   	# Out degree constraint: Total sum of x_a, with out leaving a node, is equal to 1.
   	for u in 1:T.NNodes
	    expr=0 # Start an empty expression
		for a in OutArcs[u]
			expr += x[a]
	    end
      	@constraint(TSP_model, expr==1)
    end
	
	# Constraint 5
	for a in 1:T.NArcs
		expr = 0 # Start an empty expression
		for r in AssoArcs[a]
			expr += y[r]
		end
		@constraint(TSP_model, z[a] + expr - x[a] == 0)
	end

	# Constraint 6,7,8,9
	for r in 1:T.NTriggers
		@constraint(TSP_model, y[r] - x[T.Trigger[r].trigger_arc_id] <= 0)
		@constraint(TSP_model, u[T.Arc[T.Trigger[r].trigger_arc_id].u] - u[T.Arc[T.Trigger[r].target_arc_id].u] + M*y[r] <= -1 + M)
		@constraint(TSP_model, u[T.Arc[T.Trigger[r].target_arc_id].u] - u[T.Arc[T.Trigger[r].trigger_arc_id].u] + M*w[r] <= -1 + M)
		@constraint(TSP_model, x[T.Trigger[r].trigger_arc_id] + x[T.Trigger[r].target_arc_id] + z[T.Trigger[r].target_arc_id] - w[r] <= 2)
	end

	# Constraint 10
	for a in 1:T.NArcs
		for r1 in AssoArcs[a]
			for r2 in AssoArcs[a]
				if r1 == r2
					continue
				end
				a1 = u[T.Arc[T.Trigger[r2].trigger_arc_id].u]
				a2 = M*w[r2]
				a3 = u[T.Arc[T.Trigger[r1].trigger_arc_id].u]
				a4 = y[r1]
				a5 = x[T.Trigger[r2].trigger_arc_id]
				@constraint(TSP_model, a1 - a2 -a3 + M*a4 + M*a5 <= 2*M - 1)
			end
		end
	end

   	@objective(TSP_model, Min,sum(T.Trigger[r].cost * y[r] for r in 1:T.NTriggers) + sum(T.Arc[a].cost * z[a] for a in 1:T.NArcs))

	#print(TSP_model)

   	optimize!(TSP_model)
        if (termination_status(TSP_model) != MOI.OPTIMAL)
            println(termination_status(TSP_model))
            println("It was not possible to obtain optimal solution.")
			return
        end

	if (verbose_mode)
	   println("Obtained Optimal Solution")
   	   primal_status(TSP_model)
	   
	   # --------------------------------------------------
	   if (T.NNodes > maximum_number_of_nodes_to_print)
	      println("Verbose mode: Formulation is printed only for graphs up to ",maximum_number_of_nodes_to_print," nodes.")
	   else
	      PrintHeader("Formulacao MTZ para o TSP")
   	      # print(TSP_model)
   	      # println("\nO modelo foi gravado no arquivo TSP_model.lp")
   	      # write_to_file(TSP_model, "TSP_model.lp")
	   end

	   
	   # --------------------------------------------------
	   PrintHeader("Statistics")
	   
	   ExecutionTime = time() - StartingTime
	   println("Tempo maximo de execucao: ",MaxTime," seconds.")
	   println("Tempo de execucao usado:  ",ExecutionTime," seconds.")
	   
	   NumberOfNodes = MOI.get(TSP_model, MOI.NodeCount())
	   println("Number of nodes in the Branch and Bound tree: ",NumberOfNodes)
	   
	   obj_val = objective_value(TSP_model)
	   println("Value of the solution: ",obj_val)

	   bound_val = objective_bound(TSP_model)
	   println("Value of the bound: ", bound_val)

	   ub_mtz::Vector{IntType} = Vector{IntType}(undef,T.NNodes)
	   for u in 1:T.NNodes
      	       ub_mtz[u] = -1  # indice do arco da solucao que sai do vertice u
   	   end
	   index = 1
   	   for a in 1:T.NArcs
			x_a = value(x[a])
			#println("x_", a, ": ", value(x[a]))
			if (x_a>0.0001)
				ub_mtz[index] = a
				index = index + 1
			end
	   end

	   # --------------------------------------------------
	   if (T.NNodes > maximum_number_of_nodes_to_print)
	      println("Verbose mode: Solution is printed only for graphs up to ",maximum_number_of_nodes_to_print," nodes.")
	   else
		#=
		PrintHeader("Solution of the model (variables):")
   	   	for r in 1:T.NTriggers
			println("y_",r," = ",value(y[r]))
   	   	end
		for r in 1:T.NTriggers
			println("w_",r," = ",value(w[r]))
   	   	end
		for a in 1:T.NArcs
			z_a = value(z[a])
			if (z_a!=0)
				println("z_",T.Arc[a].u,"_",T.Arc[a].v," = ",z_a)
			end
	   	end
		=#
	   end
	   
	   # --------------------------------------------------
	   PrintHeader("Solution as a sequence of nodes of the TSP cycle")
	   #=
	   u1 = 1
	   for i in 1:T.NNodes
	       print(u1,", ")
	       u1 = T.Arc[ub_mtz[u1]].v
	   end
	   println()
	   =#

	   # --------------------------------------------------
	   PrintHeader("Solution as a sequence of arcs of the TSP cycle (arc_id,u,v)")
	   for i in 1:T.NNodes
	       #println("(",ub_mtz[i],",",T.Arc[ub_mtz[i]].u,",",T.Arc[ub_mtz[i]].v,") - u[", i, "] = ", value(u[i]))
	   end
	   println()

	   # --------------------------------------------------
	   PrintHeader("Solution as arcs")
	   for a in 1:T.NArcs
			x_a = value(x[a])
			if (x_a!=0)
				println("x_",T.Arc[a].u,"_",T.Arc[a].v," = ",x_a)
			end
	   	end
	   println()
	   #println("Node 1 is in a cycle of size ",cyclesize)
	   PrintHeader("")
	end
end
