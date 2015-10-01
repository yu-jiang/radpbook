% main
filledCircle([0,-0.1],0.0075,1000,'r');
hold on
filledCircle([0,0.1],0.0075,1000,'r');
filledCircle([-0.1,0],0.0075,1000,'r');
filledCircle([0.1,0],0.0075,1000,'r');
filledCircle([0.1,.1]*sqrt(2)/2,0.0075,1000,'r');
filledCircle([0.1,-.1]*sqrt(2)/2,0.0075,1000,'r');
filledCircle([-0.1,.1]*sqrt(2)/2,0.0075,1000,'r');
filledCircle(-[0.1,.1]*sqrt(2)/2,0.0075,1000,'r');
axis([-1 1 -1 1]*.175);