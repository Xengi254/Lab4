function expmat = ECE569_MatrixLog6(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a transformation matrix T in SE(3).
% Returns the corresponding se(3) representation of exponential 
% coordinates.
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% expmat = MatrixLog6(T)
% 
% Output:
% expc6 =
%         0         0         0         0
%         0         0   -1.5708    2.3562
%         0    1.5708         0    2.3562
%         0         0         0         0

[R, p] = ECE569_TransToRp(T);
omgmat = ECE569_MatrixLog3(R);
if isequal(omgmat, zeros(3))
    expmat = [0,0,0,p(1);0,0,0,p(2);0,0,0,p(3);0,0,0,0];
else
    theta = acos((trace(R) - 1)/2);

    omg = [omgmat(3,2); omgmat(1,3); omgmat(2,1)];
    omg = omg / norm(omg);

    omg_hat = [   0       -omg(3)   omg(2)
              omg(3)     0       -omg(1)
             -omg(2)  omg(1)      0    ];

    A=(([1,0,0;0,1,0;0,0,1]-R)*omg_hat+theta*omg*transpose(omg));
    v=A\p;

    v=v*theta;
    expmat = [omgmat,v;0,0,0,0];
end
end