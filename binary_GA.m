clear;

obj_func = input('Enter the Objective_function in terms of x1 and x2 : ','s');                %taking objective function from user
f = inline(obj_func,'x1','x2');

%specifying ranges of decision variable
xmin = 0; xmax = 0.5;
ymin = 0; ymax = 0.5;

N = input('Enter the population size for solution');                    %taking population size from user
pc = input('Enter the crossover probability : ');                       %taking crossover probability from user
pm = input('Enter the mutation probability : ');                        %taking mutation probability from user

gen_max = 200;                                              %specifying maximum number of generation
str_len_x = 20;                                             %mentioning string length for x1 decision variable
str_len_y = 20;                                             %mentioning string length for x2 decision variable
str_len = str_len_x + str_len_y;

%initializing arrays required in the algorithm with all elements as zero

avg_fitness = zeros(gen_max,1);
min_fitness = zeros(gen_max,1);
max_fitness = zeros(gen_max,1);
avg_optimum_x = zeros(gen_max,1);
avg_optimum_y = zeros(gen_max,1);

iter = 0; 
pop_arr = randi([0,1],N,str_len);                           %generating initial population of size N randomly

%running while loop upto maximum no. of generation
while iter<gen_max
    
    iter = iter + 1;
    %initializing arrays with zero
    decoded_x = zeros(N,1);
    decoded_y = zeros(N,1);
    x_real = zeros(N,1);
    y_real = zeros(N,1);
    fitness_value = zeros(N,1);
    
    %algorithm for calculating fitness value
    for i = 1:N
        for j = 1:str_len_x
            decoded_x(i,1) = decoded_x(i,1) + pop_arr(i,j)*2^(str_len_x - j);
            decoded_y(i,1) = decoded_y(i,1) + pop_arr(i,j+str_len_x)*2^(str_len_x - j);
        end    
        x_real(i) = xmin + ((xmax - xmin)/(2^str_len_x - 1)) * decoded_x(i,1);
        y_real(i) = ymin + ((ymax - ymin)/(2^str_len_x - 1)) * decoded_y(i,1);
        fitness_value(i) = 1/(1+f(x_real(i),y_real(i)));
    end
    
    %roulette wheel selection(Reproduction)
    s = sum(fitness_value);
    min_fitness(iter) = min(fitness_value);
    max_fitness(iter) = max(fitness_value);
    avg_fitness(iter) = s/N;
    fitness_probability = zeros(N,1);
    cum_probability = zeros(N,1);
    mating_pool = zeros(N,str_len);

    fitness_probability(1) = fitness_value(1)/s;
    cum_probability(1) = fitness_probability(1);
    for i = 2:N
        fitness_probability(i) = fitness_value(i)/s;
        cum_probability(i) = fitness_probability(i) + cum_probability(i-1);    
    end
    
    
    
    for i = 1:N
            r = rand;
            if r < cum_probability(1)
                mating_pool(i,:) = pop_arr(1,:);
            else 
                for j = 2 :N
                    if (r <= cum_probability(j)) && (r > cum_probability(j-1))
                        mating_pool(i,:) = pop_arr(j,:);
                    end
                end
            end
    end
    
    %Crossover(two point crossover)
    mating_pair = zeros(N/2,str_len*2);
    mating_pair_1 = zeros(N/2,str_len);
    mating_pair_2 = zeros(N/2,str_len);

    random_for_pairs = randperm(N,N);
    r1 = random_for_pairs(1:N/2);
    r2 = random_for_pairs(N/2:end);
    mating_child_1 = zeros(N/2,str_len);
    mating_child_2 = zeros(N/2,str_len);
    
    for i = 1:N/2
         mating_pair(i,:) = horzcat(mating_pool(r1(i),:),mating_pool(r2(i),:)); 
         mating_pair_1(i,:) = mating_pair(i,1:str_len);
         mating_pair_2(i,:) = mating_pair(i,str_len+1:end);
         
         rand_pc = rand;
             if rand_pc <= pc
                 rand_cross = randi([1,str_len-1],2,1);
                 v = mating_pair_1(i,rand_cross(1):rand_cross(2));
                 mating_pair_1(i,rand_cross(1):rand_cross(2)) = mating_pair_2(i,rand_cross(1):rand_cross(2));
                 mating_pair_2(i,rand_cross(1):rand_cross(2)) = v;    
             end
    end
    
    mating_child = vertcat(mating_pair_1,mating_pair_2);
       
    % bitwise mutation
    for i = 1:N
        for j = 1:str_len
            rand_mut = rand;
            if rand_mut <= pm
                if mating_child(i,j) == 0
                    mating_child(i,j) = 1;
                else 
                    mating_child(i,j) = 0;
                end
            end
        end    
    end
   pop_arr = mating_child;
   optimum_x = zeros(N,1);
   optimum_y = zeros(N,1);
   x_opt = zeros(N,1);
   y_opt = zeros(N,1);
   
   %again finding real values of decision variable after performing GA
   %algorithm
   
   for i = 1:N
        for j = 1:str_len_x
            optimum_x(i,1) = optimum_x(i,1) + pop_arr(i,j)*2^(str_len_x - j);
            optimum_y(i,1) = optimum_y(i,1) + pop_arr(i,j+str_len_x)*2^(str_len_x - j);
        end  
        x_opt(i) = xmin + ((xmax - xmin)/(2^str_len_x - 1)) * optimum_x(i,1);
        y_opt(i) = ymin + ((ymax - ymin)/(2^str_len_x - 1)) * optimum_y(i,1);
   end    
   avg_optimum_x(iter) = sum(x_opt)/N;
   avg_optimum_y(iter) = sum(y_opt)/N;
   
   mating_child = zeros(N,str_len);
end

%plotting of graphs
iteration = 1:200;
figure(1)
    plot(iteration',min_fitness);
    hold on
    plot(iteration',avg_fitness);
    hold on
    plot(iteration',max_fitness);    
    hold off
    xlim([0 200]);
    ylim([0.5 1]);
    xlabel(' No. of generations ');
    ylabel(' fitness value ');
    title('fitness value vs no. of generations');
    legend({'minimum fitness', 'average fitness', 'maximum fitness'}, 'Location','southeast');
    

figure(2)
    plot(iteration',avg_optimum_x);
    hold on
    plot(iteration',avg_optimum_y);
    hold off
    xlim([0 200]);
    ylim([0 0.5]);
    xlabel(' No. of generations ');
    ylabel(' optimal solutions ');
    title('optimal solution vs no. of generations');
    legend('x1', 'x2');

    

