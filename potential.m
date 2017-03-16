% potential function accepts one or N coordinates
function V = potential(x_coord)

	if (nargin != 1)
		print_usage ("potential (x-coord)");
	end

	C(1)=-0.466516;
	C(2)=-0.834376;
	C(3)=-0.714529;
	C(4)=-0.0245586;
	C(5)=0.238837;
	C(6)=0.0143649;
	C(7)=0.271003;
	C(8)=-0.374538;
	C(9)=0.873564;
	C(10)=-0.370258;
	C(11)=0.891462;
	C(12)=-0.665239;
	C(13)=0.810546;
	C(14)=0.198216;
	C(15)=-0.816637;
	C(16)=-0.195351;
	C(17)=-0.573181;
	C(18)=0.251745;
	C(19)=0.647615;
	C(20)=0.201654;

	L=10;
	V=zeros(size(x_coord));
	for i=1:20
		V=V+C(i)*sin(2*pi*i*x_coord/L);
	end

	V(x_coord < 0) = 0;
	V(x_coord > L) = 0;
end
