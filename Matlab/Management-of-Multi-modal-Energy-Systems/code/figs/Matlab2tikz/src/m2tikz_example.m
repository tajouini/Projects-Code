figure(1),clf
plot(rand(10),rand(10))
xlabel('$x$')
ylabel('$y$')
legend('random')
grid on
box on

matlab2tikz('test.tex', 'parseStrings',false)