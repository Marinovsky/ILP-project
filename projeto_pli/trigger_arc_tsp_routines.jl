#
# Este deve ser o unico arquivo que deve ser atualizado pelo aluno
#
# Abaixo voce encontra as rotinas vazias das seguintes funcoes:
#                 TriggerArcTSP_lb_lp(T)
#                 TriggerArcTSP_lb_rlxlag(T)
#                 TriggerArcTSP_lb_colgen(T) - Opcional
#                 TriggerArcTSP_ub_lp(T)
#                 TriggerArcTSP_ub_rlxlag(T)
#                 TriggerArcTSP_ub_colgen(T) - Opcional
#                 TriggerArcTSP_ilp(T)
#
using Dates
using JuMP
using Gurobi

# --------------------------------------------------------------
function TriggerArcTSP_lb_lp(T::TriggerArcTSP, MaxTime::IntType)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here
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
   	@variable(TSP_model, 0 <= x[i=1:T.NArcs] <= 1)
	for a in 1:T.NArcs
	    set_name(x[a], "x_"*string(T.Arc[a].u)*"_"*string(T.Arc[a].v) )
	end

	# y_r
	@variable(TSP_model, 0 <= y[i=1:T.NTriggers] <= 1)
	for r in 1:T.NTriggers
	    set_name(y[r], "y_"*string(T.Trigger[r].trigger_arc_id)*"_"*string(T.Trigger[r].target_arc_id) )
	end

	# z_a
	@variable(TSP_model, 0 <= z[i=1:T.NArcs] <= 1)
	for a in 1:T.NArcs
		set_name(z[a], "z_"*string(T.Arc[a].u)*"_"*string(T.Arc[a].v) )
	end

	# w_r
	@variable(TSP_model, 0 <= w[i=1:T.NTriggers] <= 1)
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
	end
	T.time_lb_lp = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_lb_rlxlag(T::TriggerArcTSP, MaxTime::IntType)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here
	MAX_ITER = 100
	theta::Float64 = 0.3
	LB = Array{Float64}(undef, 0)

	Z_LB = 0
	Z_lambda = Inf

	sol_nodes::Vector{IntType} = Vector{IntType}(undef, T.NArcs)
	for u in 1:T.NNodes
		sol_nodes[u] = -1
	end

	k_best = 0
	M = T.NNodes + 1
	lambda::Float64 = 0.000000001

	for k=1:MAX_ITER
		(InArcs,OutArcs) = Build_In_Out_Arcs(T)
		AssoArcs = Build_Associated_Arcs(T)
		
		TSP_model = Model(Gurobi.Optimizer)
		set_optimizer_attribute(TSP_model, "OutputFlag", 1)
		
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
#=	
		# Total sum of x_a is equal to the number of nodes
		expr = 0 # Start an empty expression
		for a in 1:T.NArcs
			expr += x[a]
		end
		@constraint(TSP_model, expr==T.NNodes)
	
		# Constraint 2
		for a in 1:T.NArcs
			if T.Arc[a].v == 1
				continue
			end
			i=T.Arc[a].u
			j=T.Arc[a].v
			@constraint(TSP_model, u[i] - u[j] + M*x[a]<= -1 + M)
		end
=#
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
			#@constraint(TSP_model, y[r] - x[T.Trigger[r].trigger_arc_id] <= 0)
			#@constraint(TSP_model, u[T.Arc[T.Trigger[r].trigger_arc_id].u] - u[T.Arc[T.Trigger[r].target_arc_id].u] + M*y[r] <= -1 + M)
			#@constraint(TSP_model, u[T.Arc[T.Trigger[r].target_arc_id].u] - u[T.Arc[T.Trigger[r].trigger_arc_id].u] + M*w[r] <= -1 + M)
			#@constraint(TSP_model, x[T.Trigger[r].trigger_arc_id] + x[T.Trigger[r].target_arc_id] + z[T.Trigger[r].target_arc_id] - w[r] <= 2)
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

		# Dualize restriction 2
		lg_expr1 = 0
		for a in 1:T.NArcs
			if T.Arc[a].v != 1
				lg_expr1 += u[T.Arc[a].u] + 1 - u[T.Arc[a].v] - M * (1 - x[a])
			end
		end

		# Dualize restriction 1
		lg_expr1 = lg_expr1 + (sum(x[a] for a in 1:T.NArcs) - T.NNodes)

		# Duealize restriction 7,8,9 and 10
		for r in 1:T.NTriggers
			lg_expr1 = lg_expr1 + (y[r] - x[T.Trigger[r].trigger_arc_id])
			lg_expr1 = lg_expr1 + (u[T.Arc[T.Trigger[r].trigger_arc_id].u] - u[T.Arc[T.Trigger[r].target_arc_id].u] + M*y[r] + 1 - M)
			lg_expr1 = lg_expr1 + (u[T.Arc[T.Trigger[r].target_arc_id].u] - u[T.Arc[T.Trigger[r].trigger_arc_id].u] + M*w[r] + 1 - M)
			lg_expr1 = lg_expr1 + (x[T.Trigger[r].trigger_arc_id] + x[T.Trigger[r].target_arc_id] + z[T.Trigger[r].target_arc_id] - w[r] - 2)
		end

		# Dualize restriction 7
		#lg_expr1 = lg_expr1 + (sum(y[r] - x[T.Trigger[r].trigger_arc_id] for r in 1:T.NTriggers))

		# Dualize restriction 8
		#lg_expr1 = lg_expr1 + (sum(u[T.Arc[T.Trigger[r].trigger_arc_id].u] - u[T.Arc[T.Trigger[r].target_arc_id].u] + M*y[r] + 1 - M for r in 1:T.NTriggers))

		# Dualize restriction 9
		#lg_expr1 = lg_expr1 + (sum(u[T.Arc[T.Trigger[r].target_arc_id].u] - u[T.Arc[T.Trigger[r].trigger_arc_id].u] + M*w[r] + 1 - M for r in 1:T.NTriggers))

		# Dualize restriction 10
		#lg_expr1 = lg_expr1 + (sum(x[T.Trigger[r].trigger_arc_id] + x[T.Trigger[r].target_arc_id] + z[T.Trigger[r].target_arc_id] - w[r] - 2 for r in 1:T.NTriggers))
		try
			@objective(TSP_model, Min,sum(T.Trigger[r].cost * y[r] for r in 1:T.NTriggers) + sum(T.Arc[a].cost * z[a] for a in 1:T.NArcs) + lambda * lg_expr1)

			optimize!(TSP_model)
			if (termination_status(TSP_model) != MOI.OPTIMAL)
				println(termination_status(TSP_model))
				print("It was not possible to obtain optimal solution in iteration ", k, ".")
			end
		catch e
			break
		end
		new_Z_lambda = objective_value(TSP_model)
		
		if abs(new_Z_lambda - Z_lambda) < 0.001
			break
		else
			Z_lambda = new_Z_lambda
		end 
		append!(LB, Z_lambda)

		# Compute t for the restriction 2
		lg_expr2 = 0
		lg_expr3 = 0
		for a in 1:T.NArcs
			if T.Arc[a].v != 1
				lg_expr2 += (value(u[T.Arc[a].u]) + 1 - value(u[T.Arc[a].v]) - M * (1 - value(x[a])))^2
				lg_expr3 += (value(u[T.Arc[a].u]) + 1 - value(u[T.Arc[a].v]) - M * (1 - value(x[a])))
			end
		end

		# Compute t for the restriction 1
		lg_expr2 = lg_expr2 + (sum(value(x[a]) for a in 1:T.NArcs) - T.NNodes)^2
		lg_expr3 = lg_expr3 + (sum(value(x[a]) for a in 1:T.NArcs) - T.NNodes)

		# Compute t for the restriction 7,8,9 and 10
		for r in 1:T.NTriggers
			lg_expr2 = lg_expr2 + (value(y[r]) - value(x[T.Trigger[r].trigger_arc_id]))^2
			lg_expr3 = lg_expr3 + (value(y[r]) - value(x[T.Trigger[r].trigger_arc_id]))

			lg_expr2 = lg_expr2 + (value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) - value(u[T.Arc[T.Trigger[r].target_arc_id].u]) + M*value(y[r]) + 1 - M)^2
			lg_expr3 = lg_expr3 + (value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) - value(u[T.Arc[T.Trigger[r].target_arc_id].u]) + M*value(y[r]) + 1 - M)

			lg_expr2 = lg_expr2 + (value(u[T.Arc[T.Trigger[r].target_arc_id].u]) - value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) + M*value(w[r]) + 1 - M)^2
			lg_expr3 = lg_expr3 + (value(u[T.Arc[T.Trigger[r].target_arc_id].u]) - value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) + M*value(w[r]) + 1 - M)

			lg_expr2 = lg_expr2 + (value(x[T.Trigger[r].trigger_arc_id]) + value(x[T.Trigger[r].target_arc_id]) + value(z[T.Trigger[r].target_arc_id]) - value(w[r]) - 2 )^2
			lg_expr3 = lg_expr3 + (value(x[T.Trigger[r].trigger_arc_id]) + value(x[T.Trigger[r].target_arc_id]) + value(z[T.Trigger[r].target_arc_id]) - value(w[r]) - 2)
		end

		# Compute t for the restriction 7
		#lg_expr2 = lg_expr2 + (sum(value(y[r]) - value(x[T.Trigger[r].trigger_arc_id]) for r in 1:T.NTriggers))^2
		#lg_expr3 = lg_expr3 + (sum(value(y[r]) - value(x[T.Trigger[r].trigger_arc_id]) for r in 1:T.NTriggers))

		# Compute t for the restriction 8
		#lg_expr2 = lg_expr2 + (sum(value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) - value(u[T.Arc[T.Trigger[r].target_arc_id].u]) + M*y[r] + 1 - M for r in 1:T.NTriggers))^2
		#lg_expr3 = lg_expr3 + (sum(value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) - value(u[T.Arc[T.Trigger[r].target_arc_id].u]) + M*y[r] + 1 - M for r in 1:T.NTriggers))

		# Compute t for the restriction 9
		#lg_expr2 = lg_expr2 + (sum(value(u[T.Arc[T.Trigger[r].target_arc_id].u]) - value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) + M*w[r] + 1 - M for r in 1:T.NTriggers))^2
		#lg_expr3 = lg_expr3 + (sum(value(u[T.Arc[T.Trigger[r].target_arc_id].u]) - value(u[T.Arc[T.Trigger[r].trigger_arc_id].u]) + M*w[r] + 1 - M for r in 1:T.NTriggers))

		# Compute t for the restriction 10
		#lg_expr2 = lg_expr2 + (sum(value(x[T.Trigger[r].trigger_arc_id]) + value(x[T.Trigger[r].target_arc_id]) + value(z[T.Trigger[r].target_arc_id]) - value(w[r]) - 2 for r in 1:T.NTriggers))^2
		#lg_expr3 = lg_expr3 + (sum(value(x[T.Trigger[r].trigger_arc_id]) + value(x[T.Trigger[r].target_arc_id]) + value(z[T.Trigger[r].target_arc_id]) - value(w[r]) - 2 for r in 1:T.NTriggers))
		
		t = (theta^(k))*(Z_LB - Z_lambda) / lg_expr2

		println()
		println("theta ", theta)
		println("Z_LB ", Z_LB)
		println("Z_lambda ", Z_lambda)
		println("t_ ", k, ": ", t)
		lambda = max(lambda + t*(lg_expr3), 0)
		println("lambda ", lambda)

		if max(Z_LB, Z_lambda) == Z_lambda
			Z_LB = Z_lambda
			k_best = k
			for a in 1:T.NArcs
				x_a = value(x[a])
				if (x_a>0.0001)
				  sol_nodes[a] = x_a
				end
			end
		end
		empty!(x)
		empty!(y)
		empty!(z)
		empty!(w)
	end
	ExecutionTime = time() - StartingTime # To compute the time used by this routine.
	println("Tempo de execucao usado:  ", ExecutionTime," seconds.")
	println("Best value of the solution: ", Z_LB)
	println("In iteration: ", k_best)
	PrintHeader("Solution as a sequence of nodes of the TSP cycle")
	PrintHeader("Solution as arcs")
	for a in 1:T.NArcs
		x_a = value(sol_nodes[a])
		if (x_a!=0)
			#println("x_",T.Arc[a].u,"_",T.Arc[a].v," = ",x_a)
		end
	end
	PrintHeader("List of lower bounds")
	for i in 1:length(LB)
		println("[",i,"]: ", LB[i])
	end
end
# --------------------------------------------------------------
function TriggerArcTSP_lb_colgen(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_lb_colgen = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_lp(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ub_lp = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_rlxlag(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ub_rlxlag = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_colgen(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ub_colgen = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ilp(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ilp = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
