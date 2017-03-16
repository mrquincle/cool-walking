clf

L=10;
x=0:0.001:L;

V=zeros(size(x));

V=potential(x);

plot(x,V, '-b', 'LineWidth',4);
hold on

k=1.38064852*10^-23;
k=1; % apparently set to 1 in this case: beta = 1/(kT)

% high temperature
T=3;
beta=1/(k*T);
%output_precision(32)

rho_un=exp(-beta * V);
Z=sum(rho_un)/size(x,2)*L;
Z

rho=rho_un/Z;

plot(x,rho, '-r', 'LineWidth',4);

% low temperature
T=0.1;
beta=1/(k*T);

rho_un=exp(-beta * V);
max(rho_un)
Z=sum(rho_un)/size(x,2)*L;
Z

rho=rho_un/Z;
max(rho)

plot(x,rho,'-g', 'LineWidth',4);

title("One-dimensional potential energy V(x)")
set(gca, "linewidth", 4, "fontsize", 30);
L = legend("V(x)", "T = 3.0", " T = 0.1");
legend("boxoff");
xlabel("x");
ylabel("V(x) or \\rho(x)")
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize', 12);

set(L,'FontSize', 8);
hold off;

print(h,'-dpng','-color','potential_energy.png')

