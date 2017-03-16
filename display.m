clear all
clf

% structured
load distances.mat
c1 = convergence;

% uninformed
load distances3.mat;
c2 = convergence;

% J-walking
load distances2.mat;
c3 = convergence;

h = figure(1);

plot(c1,'-g','LineWidth',4);
hold on;
plot(c2,'-b','LineWidth',4);
plot(c3,'-r','LineWidth',4);

title("Compare J-Walking with Hybrid Proposals")
set(gca, "linewidth", 4, "fontsize", 30);
L = legend("Structured hybrid proposals", "Uninformed hybrid proposals", "J-Walking");
legend("boxoff");
xlabel("t (MH steps x 1000)");
ylabel("\\chi(t)")
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize', 12);

set(L,'FontSize', 8);
hold off;

print(h,'-dpng','-color','mcmc_hybrids.png')

