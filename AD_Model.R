AD_Model = function(n=500, #number of agents
                    k=6, #number of preference dimensions
                    d=.75, #decay parameter (lower D = faster decay; vary in range 0-1)
                    iters = 5000000, #number of ticks the model will run for
                    measurement_interval = 100000, #how many ticks between population measurements
                    SmallWorld = FALSE, #T = small-world topology, F = fully-connected net
                    caves = 1, #number of caves in small-world; doesn't affect fully-connected nets
                    rewire_proportion = .10) #proportion of rewired edges in small-world; doesn't affect fully-connected nets
{
  
  #Assign agent attributes:
  
  #1 - Preferences
  #Create data frame with one row per agent and one column/variable per preference
  #Each cell will be given a random value between -1 and 1  
  p = data.frame(matrix(runif(n*k, min= -1, max=1), nrow = n, ncol = k))
  
  #2 - Associations
  #Create an empty list
  a = list()
  
  #For each agent, we create a K x K matrix (K = number of preference dimensions) 
  #Each cell initially has the value 1 in it
  #Add each matrix to the list 
  mat = matrix(1, nrow = k, ncol = k)
  for(y in 1:n) {
    a[[y]] = mat
  }
  
  #Create data frame to store results
  tick = seq(1:iters)
  results = data.frame(tick) #store results here
  results$pref_sim=NA #preference similarity (see Goldberg & Stein for explanations of measures)
  results$pref_cong=NA #preference congruence 
  
  #Create network topology
  if(SmallWorld == T) {
    cavesize = n / caves
    net = make_empty_graph(directed=T)
    for(c in 1:caves) {
      net2 = make_full_graph(cavesize, directed=T, loops=F)
      V(net2)$cave = c
      net = net + net2
    }
    
    #rewire network
    net = rewire(net, with = keeping_degseq(niter = (length(E(net))*rewire_proportion)) ) #rewire edges
    net = simplify(net)
    net = delete_vertices(net, degree(net)==0)
  }
  
  #Update agents
  for(i in 1:iters) {
    
    #1 - An actor and observer (so 2 agents total) are randomly selected to interact
    if(SmallWorld == T) {
      #small world topology 
      observer = sample(1:n, 1) 
      actor = sample(adjacent_vertices(net, observer, mode=c("out"))[[1]], 1)
    } else{ #Fully connected
      agents = sample(1:n, 2) #First is actor, Second is observer 
      actor = agents[1]
      observer = agents[2]
    }
    
    #2- The actor (Agent 1) enacts two practices (based on her preferences)
    behaviors = sample(1:k, #1 through K 
                       2, #pick two preferences
                       prob = (exp( p[actor, ])/( sum(exp(p[actor, ])))))
    
    #3 - The observer (Agent 2) updates perception of associations between these practices
    a[[observer]][behaviors[1], behaviors[2]] = a[[observer]][behaviors[1], behaviors[2]] + 1
    a[[observer]][behaviors[2], behaviors[1]] = a[[observer]][behaviors[2], behaviors[1]] + 1 
    
    #4 - Observer (Agent 2) updates one preference (the weaker of the two)
    #"Weaker" here means that the absolute value of the distance between observer's preference and the mean preference is lower
    updated_preference = ifelse((abs(p[observer, behaviors[1]] - rowMeans(p[observer,]))) < (abs(p[observer, behaviors[2]] - rowMeans(p[observer,]))), 
                                behaviors[1], #return weaker of two preferenceds
                                behaviors[2])
    
    #5 - Observer updates selected preference from step 4 
    old_prefs = p[observer,]
    
    new_prefs = old_prefs
    #new_prefs[updated_preference] = old_prefs[updated_preference] + runif(1, min = -1, max = 1)
    #update using a normal instead of uniform distribution
    new_prefs[updated_preference] = old_prefs[updated_preference] + rnorm(1, mean = 0, sd = 1)
    
    #6 - Test constraint satisfaction
    
    #convert observer's preference vector into matrix of differences
    old_pref_matrix = abs(outer(as.numeric(old_prefs), as.numeric(old_prefs), FUN="-"))
    new_pref_matrix = abs(outer(as.numeric(new_prefs), as.numeric(new_prefs), FUN="-"))
    
    #standarize all three relevant matrices (old pref, new pref, and associations) by maximum values
    old_pref_matrix = old_pref_matrix/max(old_pref_matrix)
    new_pref_matrix = new_pref_matrix/max(new_pref_matrix)
    
    assoc_std = a[[observer]]/max(a[[observer]])
    
    #Calculate CS for old and new preferences
    old_CS = (1 / (k * (k - 1))) * sum(abs(assoc_std - old_pref_matrix))
    new_CS = (1 / (k * (k - 1))) * sum(abs(assoc_std - new_pref_matrix)) 
    
    #If CS with new preference is higher: Keep NEW preference
    #Otherwise: Keep OLD preference
    if(new_CS > old_CS) {
      p[observer,]=new_prefs
    }
    
    #7 - Apply decay function to associations
    a[[observer]] = a[[observer]] * d
    
    #Record population-level outcomes
    if( (i/measurement_interval) == round(i/measurement_interval) | i == 1) {
      print(i)
      pref_sim_mat = matrix(NA, nrow=n, ncol=n)
      for(x in 1:n) {
        for(y in 1:n) {
          pref_sim_mat[x, y]=cor(as.numeric(p[x,]), as.numeric(p[y,]))
        }
      }
      pref_sim_mat[lower.tri(pref_sim_mat)] = NA #preference similarity only measured for unique pairs
      diag(pref_sim_mat)=NA #set diagonals to not count
      results$pref_sim[i] = mean(pref_sim_mat, na.rm=T)
      results$pref_cong[i] = mean(abs(pref_sim_mat), na.rm=T)  
      
      #Stop if model has converged
      if(results$pref_cong[i] >= .99) {
        break
      }
    }
  }
  
  #Clean and return results data frame
  results = subset(results, !(is.na(results$pref_cong)) )
  return(results)
}
