function [d dp] = displacements()
	x=0:0.001:10;
	V = -potential(x);
	V(V<3)=0;
	
	[peaks locs] = findpeaks(V);
	dp = x(locs);
	
	%plot(x,V,'-r')
	%hold on
	%plot(dp,zeros(size(dp)),'o');
	%hold off;

	% form matrix, extract upper triangular as packed column and take absolute
	d = abs(triu(repmat(dp,size(dp'))-repmat(dp',size(dp)),1,"pack"));
	d = [d; -d];
end

