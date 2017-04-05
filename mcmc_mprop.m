% MCMC 

% use multiple proposals, rather than temperature

clear all;
clf;

% steps (takes a minute or two for 200 000 steps)
T = 200000;

% burning
B = 1000;

% thinning
S = 10;

% convergence
R = 1000;

convergence = zeros(1,T/R);

% In this proposal we do not use a Hamiltonian, but simple local updates. It is possible to include Hamiltonian
% updates as well later. It is an independent improvement with respect to the implemented methods.
function new_state = propose(state, jump, method) 
	if (jump)
		switch(method)
			case "informed_mixture_walking"
				d = displacements();
				di = d(randi(length(d)));
				proposal = di + normrnd(0,0.1); % large jumps, with info about displacements... not 3x out of 12
			otherwise
				proposal = normrnd(0,5); % just large jumps
		end
	else 
		proposal = normrnd(0,0.1); % micro-steps
	end
	new_state = state + proposal;
end

function p = pdf(state, temp)
	beta = 1/temp;
	p = exp (-beta * potential(state));
	% normalize, only precalculated for T=0.1, and T=3.0
	if (temp == 0.1)
		p = p / 1.4066242476099965000000000000000e+15;
	elseif (temp == 3)
		p = p / 1.1820461336211799618922668742016e+01;
	else
		assert(false);
	end
end

Temp_std = 0.1;

% start at random location
state(1)=5;
p(1) = pdf(state(1),Temp_std);

% we should start with non-zero prob, or else we have /0.
assert(p(1) != 0);
		
succ_jumps = 0;
fail_jumps = 0;

chi_resolution=0.0001;
true_resolution=0.000001;

L=10;

% first reshape in 100x100 matrix and then average over rows, to be able to compare with histogram
% href / 10000*10 is to normalize it
%href=(sum(pdf(reshape([true_resolution:true_resolution:L],chi_resolution/(true_resolution),10/chi_resolution)',0.1),2)*true_resolution/chi_resolution);
href=(sum(pdf(reshape([true_resolution:true_resolution:L],chi_resolution/(true_resolution),10/chi_resolution)',0.1),2)*true_resolution);
length(href)
sum(href)
	
% This is not "exactly" j-walking, but a more lazy implementation in which we do not run two chains in parallel,
% but jump with probability p_jump to a state in the high temperature walker. This skews the results in the 
% following way: the chain will only be able to reach high-temperature states jumping from low-temperature ones.
% Consequence: This means we might miss a mode if it is too far from the other ones. Hence it is easy to see if this 
% happens and we can disregard a run like that for our comparison.

method="jwalking";
method="jwalking_simultaneously";
method="mixture_walking";
method="informed_mixture_walking";

Temp_jump = 3;

switch(method)
	case "jwalking_simultaneously"
		hot_state(1)=5;
		hot_p(1) = pdf(hot_state(1), Temp_jump);
	otherwise
end

% Let the MCMC sampler run over T time steps. Note that we do not include HMC gradient descent steps. Our chain
% should henceforth convergence slower than in [1].
for t=1:T
	% Set the jump probability to p_jump = 0.03.
	u = rand(1);
	p_jump = 0.03;
	jump = (u < p_jump);

	switch(method)
		case "jwalking"
			state(t+1) = propose(state(t), false, method);
			if (jump)
				p(t+1) = pdf(state(t+1), Temp_jump);
			else
				p(t+1) = pdf(state(t+1), Temp_std);
			end
		case "jwalking_simultaneously"
			hot_state(t+1) = propose(hot_state(t), false, method);
			hot_p(t+1) = pdf(hot_state(t+1), Temp_jump);
			if (jump)
				state(t+1) = hot_state(t+1);
				% accept with low-temperature probability
				p(t+1) = pdf(state(t+1), Temp_std);
				%p(t+1) = hot_p(t+1);
			else
				state(t+1) = propose(state(t), jump, method);
				p(t+1) = pdf(state(t+1), Temp_std);
			end
		case "mixture_walking"
			state(t+1) = propose(state(t), jump, method);
			p(t+1) = pdf(state(t+1), Temp_std);
		case "informed_mixture_walking"
			state(t+1) = propose(state(t), jump, method);
			p(t+1) = pdf(state(t+1), Temp_std);
		otherwise
			printf "Unknown method\n"
	end
	
	alpha = p(t+1)/p(t);
	u = rand(1);
	if (alpha <= u)
		p(t+1) = p(t); %reject
		state(t+1)=state(t);
		if (jump) fail_jumps++; end
	else
		if (jump) succ_jumps++; end
	end

	if (mod(t,R) == 0) 
	%	h0 = (hist(state,chi_resolution:chi_resolution:L,1)'/chi_resolution);
		h0 = hist(state,chi_resolution:chi_resolution:L,1)';
		%h0(1:1/(chi_resolution*50):1/chi_resolution)'
	
		% over 10.000 points we calculate all differences, hence, before the sqrt we divide by chi_resolution*L
		%conv = sqrt(sum((h0-href).^2)*chi_resolution*L);
		%conv = sqrt( sum((h0-href).^2)*chi_resolution*L);
		conv = sqrt( sum((h0-href).^2));
		t
		conv
		convergence(t/R) = conv;
	end
end

%figure(3);
%bar([0.1:0.1:10],href);
%figure(4);
%bar([0.1:0.1:10],conv);

% burning
state = state(B:length(state));

% thinning
state = state(1:S:length(state));

figure(1)

hist(state,100)

st_avg = mean(state);
printf("Average <x>: %f\n", st_avg);

v_avg = mean(potential(state));
printf("Average <V>: %f\n", v_avg);

printf("Successful jumps: %i\n", succ_jumps);
printf("Failing jumps: %i\n", fail_jumps);

percentage=(succ_jumps*100)/(succ_jumps + fail_jumps);

printf("Percentage: %4f%%\n", percentage);

hold on

figure(2)

t=1:length(convergence);
plot(t, convergence,'-');


figure(3)
plot(h0,'.b')
hold on
plot(href,'.r')
title("Compare h0 with href");
hold off

figure(4)
bins=1000;
hist(state*bins/10,0:1:bins-1,true,'b')
hold on
r=sum(reshape(href,length(href)/bins,bins))';
bar(r,'hist','g')
axis([0 bins*0.6]); % interesting part
title("Compare again with bar plots")
hold off


