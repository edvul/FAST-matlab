addpath(genpath(pwd));
savepath();
rehash;

try
    test = fastStart(2,'funcVal', {[0 1] [0.001 1]});
    fprintf('Path setup completed.\n');
catch
    error('Setup failed -- something is wrong');
end