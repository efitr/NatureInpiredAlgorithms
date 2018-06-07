##Input
    #Problem_size, 
    #Cells_num, 
    #N_ed, number of elimination-dispersal steps or removal of bacteria as a result of gradual or unexpected changes
    #N_re, number of reproduction steps
    #N_c, number of chemotaxis steps
    #N_s, number of swim steps for a given cel
    #Step_size, is a random direction vector with the same number of dimensions as the problem space, each value belongs between -1 and 1
    #d_attract, attraction coeficient
    #w_attract, attraction coeficient
    #h_repellant, repel coeficient
    #w_repellant, repel coeficient
    #P_ed, is the probability of a cell beign subjected to elimination and dispersal.

##Output
    #Cell_best

#Population <-- InitializePopulation(Cell_num, 
                                    #Problem_size)

#For(l = () to N_ed)
    #For(k = () to N_re)
        #For(j = () to N_c)
            #ChemotaxisAndSwim(Population, 
                            #Problem_size, 
                            #Cells_num, 
                            #N_s, 
                            #Step_size, 
                            #d_attract, 
                            #w_attract, 
                            #h_repellant, 
                            #w_repellant)

            #For(Cell belongs Population)
                #If(Cost(cell)<=Cost(Cell_best))
                    #Cell_best <-- Cell
                #End
            #End
        #End
        #SortByCellHealth(Population)

        #Selected <-- SelectByCellHealth(Population, Cells_num/2)
        #Population <-- Selected
        #Population <-- Selected
    #End
    #For(Cell belongs Population)
        #If(Rand()<=P_ed)
            #Cell <-- CreateCellAtRandomLocation()
        #End
    #End
#End

#Return (Cell_best)

#Pseudocode 

#Input
    #Population
    #Problem_size, 
    #Cells_num, 
    #N_s, number of swim steps for a given cel
    #Step_size, is a random direction vector with the same number of dimensions as the problem space, each value belongs between -1 and 1
    #d_attract, attraction coeficient
    #w_attract, attraction coeficient
    #h_repellant, repel coeficient
    #w_repellant, repel coeficient
    
#For(Cell belongs Population)
    #Cell_fitness <-- Cost(Cell)+Interaction(Cell, 
                                            #Population, 
                                            #d_attract, 
                                            #w_attract, 
                                            #h_repellant, 
                                            #w_repellant) 

    #Cell_health <-- Cell_fitness
    #Cell' <-- null
    #For(i=0 to N_s)
        #RandomStepDirection <-- CreateStep(Problem_size)
        #Cell' <-- TakeStep(RandomStepDirection, Step_size)
        #Cell'_fitness <-- Cost(Cell')+Interaction(Cell', 
                                                #Population, 
                                                #d_attract, 
                                                #w_attract, 
                                                #h_repellant, 
                                                #w_repellant)
        #If(Cell'_fitness > Cell_fitness)
            #i <-- N_s
        #Else
            #Cell <-- Cell'
            #Cell_health <-- Cell_health, Cell'_fitness
        #End
    #End
#End

#Heuristics
    #d_attract = 0.1
    #w_attract = 0.2
    #h_repellant = d_attract
    #w_repellant = 10
    #Step_size = /*Small Fraction of the search space*/ 0.1
    #Reproduction half the population with a low health metric are discarded
    # and two copies of each member from the first (high-health) half of the population are retained
    #Prob of elimination and dispersal (P_ed) is commonly set quite large, such as 0.25

    def objective_function(vector)
        return vector.inject(0.0) {|sum, x| sum +  (x ** 2.0)}
      end
      
      def random_vector(minmax)
        return Array.new(minmax.size) do |i|
          minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * rand())
        end
      end
      
      def generate_random_direction(problem_size)
        bounds = Array.new(problem_size){[-1.0,1.0]}
        return random_vector(bounds)
      end
      
      def compute_cell_interaction(cell, cells, d, w)
        sum = 0.0
        cells.each do |other|
          diff = 0.0
          cell[:vector].each_index do |i|
            diff += (cell[:vector][i] - other[:vector][i])**2.0
          end
          sum += d * Math.exp(w * diff)
        end
        return sum
      end
      
      def attract_repel(cell, cells, d_attr, w_attr, h_rep, w_rep)
        attract = compute_cell_interaction(cell, cells, -d_attr, -w_attr)
        repel = compute_cell_interaction(cell, cells, h_rep, -w_rep)
        return attract + repel
      end
      
      def evaluate(cell, cells, d_attr, w_attr, h_rep, w_rep)
        cell[:cost] = objective_function(cell[:vector])
        cell[:inter] = attract_repel(cell, cells, d_attr, w_attr, h_rep, w_rep)
        cell[:fitness] = cell[:cost] + cell[:inter]
      end
      
      def tumble_cell(search_space, cell, step_size)
        step = generate_random_direction(search_space.size)
        vector = Array.new(search_space.size)
        vector.each_index do |i|
          vector[i] = cell[:vector][i] + step_size * step[i]
          vector[i] = search_space[i][0] if vector[i] < search_space[i][0]
          vector[i] = search_space[i][1] if vector[i] > search_space[i][1]
        end
        return {:vector=>vector}
      end
      
      def chemotaxis(cells, search_space, chem_steps, swim_length, step_size,
          d_attr, w_attr, h_rep, w_rep)
        best = nil
        chem_steps.times do |j|
          moved_cells = []
          cells.each_with_index do |cell, i|
            sum_nutrients = 0.0
            evaluate(cell, cells, d_attr, w_attr, h_rep, w_rep)
            best = cell if best.nil? or cell[:cost] < best[:cost]
            sum_nutrients += cell[:fitness]
            swim_length.times do |m|
              new_cell = tumble_cell(search_space, cell, step_size)
              evaluate(new_cell, cells, d_attr, w_attr, h_rep, w_rep)
              best = cell if cell[:cost] < best[:cost]
              break if new_cell[:fitness] > cell[:fitness]
              cell = new_cell
              sum_nutrients += cell[:fitness]
            end
            cell[:sum_nutrients] = sum_nutrients
            moved_cells << cell
          end
          puts "  >> chemo=#{j}, f=#{best[:fitness]}, cost=#{best[:cost]}"
          cells = moved_cells
        end
        return [best, cells]
      end
      
      def search(search_space, pop_size, elim_disp_steps, repro_steps,
          chem_steps, swim_length, step_size, d_attr, w_attr, h_rep, w_rep,
          p_eliminate)
        cells = Array.new(pop_size) { {:vector=>random_vector(search_space)} }
        best = nil
        elim_disp_steps.times do |l|
          repro_steps.times do |k|
            c_best, cells = chemotaxis(cells, search_space, chem_steps,
              swim_length, step_size, d_attr, w_attr, h_rep, w_rep)
            best = c_best if best.nil? or c_best[:cost] < best[:cost]
            puts " > best fitness=#{best[:fitness]}, cost=#{best[:cost]}"
            cells.sort{|x,y| x[:sum_nutrients]<=>y[:sum_nutrients]}
            cells = cells.first(pop_size/2) + cells.first(pop_size/2)
          end
          cells.each do |cell|
            if rand() <= p_eliminate
              cell[:vector] = random_vector(search_space)
            end
          end
        end
        return best
      end
      
      if __FILE__ == $0
        # problem configuration
        problem_size = 2
        search_space = Array.new(problem_size) {|i| [-5, 5]}
        # algorithm configuration
        pop_size = 50
        step_size = 0.1 # Ci
        elim_disp_steps = 1 # Ned
        repro_steps = 4 # Nre
        chem_steps = 70 # Nc
        swim_length = 4 # Ns
        p_eliminate = 0.25 # Ped
        d_attr = 0.1
        w_attr = 0.2
        h_rep = d_attr
        w_rep = 10
        # execute the algorithm
        best = search(search_space, pop_size, elim_disp_steps, repro_steps,
          chem_steps, swim_length, step_size, d_attr, w_attr, h_rep, w_rep,
          p_eliminate)
        puts "done! Solution: c=#{best[:cost]}, v=#{best[:vector].inspect}"
      end