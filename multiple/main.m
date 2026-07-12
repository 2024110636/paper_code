N_trials = 5;
fitness_all = zeros(N_trials, 1);

for t = 1:N_trials
    [fitness_all(t), ~] = run_GA_full(440 + t); % 每次不同随机种子
end

fprintf('Mean: %.4f | Std: %.4f\n', mean(fitness_all), std(fitness_all));
boxplot(fitness_all); ylabel('Best Fitness');

