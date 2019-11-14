fs = 30e9;
targetPos = [20; 20];   % target position
refPos = [-10 10 0; ... % reference receiver positions [x; y]
          -10 10 10];

x = 2*randi([0 1], 100, 1)-1;
      
[y, delays] = add_delay(x, targetPos, refPos, fs);