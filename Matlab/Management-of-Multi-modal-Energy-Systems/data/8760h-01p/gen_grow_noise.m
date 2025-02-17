function  X = gen_grow_noise_DA(N, noise, real, start)
    relGrowth = linspace(0.01,0.3,N); % from 1% to 30% noise
    relGrowth = relGrowth*sum(relGrowth)/N; % average
    l = length(real)-start+1;
    X = zeros(1,l);
    for i=1:l/N
        X((i-1)*N+start:i*N+start-1) = real((i-1)*N+start:i*N+start-1) +...
        real((i-1)*N+start:i*N+start-1).*relGrowth.*noise((i-1)*N+start:i*N+start-1);
    end
end