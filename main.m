%% load data
clear; close all; clc;
load('ElecPosXYZ');
load("ElecPatch.mat")
load("Interictal.mat")
%% part a 
% Forward Matrix
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1;
[LocMat,G] = ForwardModel_3shell(Resolution, ModelParams);
% plot dipole locations
figure();
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); axis('equal');
xlabel('X', 'Interpreter','latex')
ylabel('Y', 'Interpreter','latex')
zlabel('Z', 'Interpreter','latex')
title('dipoles positions', 'Interpreter','latex')
%% part b
ElecNames = cell(1, 21);
ElecXYZ = zeros(3, 21);
for i = 1:21
    ElecNames{i} = ElecPos{1, i}.Name;
    ElecXYZ(:, i) = ElecPos{1, i}.XYZ';
end
ElecXYZ = ElecXYZ * 9.2;
figure();
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); axis('equal');
xlabel('X', 'Interpreter','latex')
ylabel('Y', 'Interpreter','latex')
zlabel('Z', 'Interpreter','latex')
title('dipoles/electrodes positions', 'Interpreter','latex');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end
n = 1289;
e_q = LocMat(:, n) / norm(LocMat(:, n));  
line([LocMat(1,n),LocMat(1,n)+e_q(1)],[LocMat(2,n),LocMat(2,n)+e_q(2)],[LocMat(3,n),LocMat(3,n)+e_q(3)],'Color','yellow','linewidth',1.5);
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled') % selected dipole
legend('dipoles', 'electrodes', 'selected dipole')
%% part d
G_select = G(:, 3*n-2:3*n);
e_q = LocMat(:, n) / norm(LocMat(:, n));
eq_superficial = e_q;
Q = e_q * Interictal(1, :);
M = G_select*Q;
disp_eeg(M, [], 256, ElecNames);
%% part e
% find the channel with maximum voltage
max_amp = -inf;
max_idx = 0;
for i = 1:21
    if max(M(i, :) > max_amp)
        max_amp = max(M(i, :));
        max_idx = i;
    end
end
% locate windows with respect to channel with maximum peak voltage
[~, locs] = findpeaks(M(max_idx, :), 'MinPeakDistance', 10, 'MinPeakHeight', 10);
mean_signal = zeros(21,1);
for i = 1:21
    for j = 1:length(locs)
        mean_signal(i) = mean_signal(i) + mean(M(i, locs(j)-3:locs(j)+3))/length(locs);
    end
end
figure();
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 36, mean_signal, 'o', 'filled');
colorbar;
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end
%% part f
Display_Potential_3D(9.2, mean_signal)
%% part g
alpha = 0.1;
% MNE
Q_MNE = G'*(G*G' + alpha*eye(21))^(-1)*mean_signal;
% WMNE
P = size(G, 2) / 3;
N = 21;
Omega = zeros(1, P);
for i = 1:P
    temp = 0;
    for j = 1:N
        g = G(j, 3*i-2:3*i);
        temp = temp + g*g';
    end
    Omega(i) = sqrt(temp);
end
W_WMNE = kron(diag(Omega), eye(3));
Q_WMNE = (W_WMNE'*W_WMNE)^(-1)*G'*(G*(W_WMNE'*W_WMNE)^(-1)*G' + alpha*eye(N))^(-1)*mean_signal;
% LORETA
d = 1;
A1 = zeros(P, P);
for alpha1 = 1:P
    for beta = 1:P
        A1(alpha1, beta) = (norm(LocMat(:, alpha1) - LocMat(:, beta)) == d)/6;
    end
end
A0 = diag(A1*ones(P, 1))^(-1)*A1;
A = kron(A0, eye(3));
B = 6/(d^2)*(A-eye(3*P));
W_LORETA = kron(diag(Omega), eye(3)) * (B' * B) * kron(diag(Omega), eye(3));
Q_LORETA = (W_LORETA'*W_LORETA)^(-1)*G'*(G*(W_LORETA'*W_LORETA)^(-1)*G' + alpha*eye(N))^(-1)*mean_signal;
% SLORETA
S_Q = G'*(G*G' + alpha*eye(N))^(-1)*G;
Q_SLORETA = zeros(3*P, 1);
for i = 1:P
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end
%% part h
Q_MNE_max = -inf;
Q_MNE_max_idx = 0;
Q_WMNE_max = -inf;
Q_WMNE_max_idx = 0;
Q_LORETA_max = -inf;
Q_LORETA_max_idx = 0;
Q_SLORETA_max = -inf;
Q_SLORETA_max_idx = 0;
for i = 1:P
    if norm(Q_MNE(3*i-2:3*i)) > Q_MNE_max
        Q_MNE_max = norm(Q_MNE(3*i-2:3*i));
        Q_MNE_max_idx = i;
    end
    if norm(Q_WMNE(3*i-2:3*i)) > Q_WMNE_max
        Q_WMNE_max = norm(Q_WMNE(3*i-2:3*i));
        Q_WMNE_max_idx = i;
    end
    if norm(Q_LORETA(3*i-2:3*i)) > Q_LORETA_max
        Q_LORETA_max = norm(Q_LORETA(3*i-2:3*i));
        Q_LORETA_max_idx = i;
    end
    if norm(Q_SLORETA(3*i-2:3*i)) > Q_SLORETA_max
        Q_SLORETA_max = norm(Q_SLORETA(3*i-2:3*i));
        Q_SLORETA_max_idx = i;
    end
end
figure();
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); axis('equal');
xlabel('X', 'Interpreter','latex')
ylabel('Y', 'Interpreter','latex')
zlabel('Z', 'Interpreter','latex')
title('dipoles/electrodes positions', 'Interpreter','latex');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end

scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled')
scatter3(LocMat(1, Q_MNE_max_idx), LocMat(2, Q_MNE_max_idx), LocMat(3, Q_MNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_WMNE_max_idx), LocMat(2, Q_WMNE_max_idx), LocMat(3, Q_WMNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_LORETA_max_idx), LocMat(2, Q_LORETA_max_idx), LocMat(3, Q_LORETA_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_SLORETA_max_idx), LocMat(2, Q_SLORETA_max_idx), LocMat(3, Q_SLORETA_max_idx), 's', 'filled')
legend('dipoles', 'electrodes', 'selected dipole', 'predict (MNE)', 'predict (WMNE)', 'predict (LORETA)', 'predict (SLORETA)')

% location error
LocErr_MNE = mse(LocMat(:, n), LocMat(:, Q_MNE_max_idx)); %#ok<NASGU> 
LocErr_WMNE = mse(LocMat(:, n), LocMat(:, Q_WMNE_max_idx)); %#ok<NASGU> 
LocErr_LORETA = mse(LocMat(:, n), LocMat(:, Q_LORETA_max_idx)); %#ok<NASGU>
LocErr_SLORETA = mse(LocMat(:, n), LocMat(:, Q_SLORETA_max_idx)); %#ok<NASGU> 
% angle error
e_q_MNE_predict = Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx) / norm(Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx));
AngErr_MNE = acos((e_q'*e_q_MNE_predict) / (norm(e_q)*norm(e_q_MNE_predict)))*180/pi; %#ok<NASGU> 
e_q_WMNE_predict = Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx) / norm(Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx));
AngErr_WMNE = acos((e_q'*e_q_WMNE_predict) / (norm(e_q)*norm(e_q_WMNE_predict)))*180/pi; %#ok<NASGU> 
e_q_LORETA_predict = Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx) / norm(Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx));
AngErr_LORETA = acos((e_q'*e_q_LORETA_predict) / (norm(e_q)*norm(e_q_LORETA_predict)))*180/pi; %#ok<NASGU> 
e_q_SLORETA_predict = Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx) / norm(Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx));
AngErr_SLORETA = acos((e_q'*e_q_SLORETA_predict) / (norm(e_q)*norm(e_q_SLORETA_predict)))*180/pi;  %#ok<NASGU> 
%% part b
ElecNames = cell(1, 21);
ElecXYZ = zeros(3, 21);
for i = 1:21
    ElecNames{i} = ElecPos{1, i}.Name;
    ElecXYZ(:, i) = ElecPos{1, i}.XYZ';
end
ElecXYZ = ElecXYZ * 9.2;
figure();
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); axis('equal');
xlabel('X', 'Interpreter','latex')
ylabel('Y', 'Interpreter','latex')
zlabel('Z', 'Interpreter','latex')
title('dipoles/electrodes positions', 'Interpreter','latex');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end
n = 853;
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled') % dipole
legend('dipoles', 'electrodes', 'selected dipole')
%% part d
G_select = G(:, 3*n-2:3*n);
e_q = LocMat(:, n) / norm(LocMat(:, n));
e_q_deep = e_q;
Q = e_q * Interictal(1, :);
M = G_select*Q;
disp_eeg(M, [], 256, ElecNames);
%% part e
% find the channel with maximum voltage
max_amp = -inf;
max_idx = 0;
for i = 1:21
    if max(M(i, :) > max_amp)
        max_amp = max(M(i, :));
        max_idx = i;
    end
end
% locate windows with respect to channel with maximum peak voltage
[~, locs] = findpeaks(M(max_idx, :), 'MinPeakDistance', 10, 'MinPeakHeight', 10);
mean_signal = zeros(21,1);
for i = 1:21
    for j = 1:length(locs)
        mean_signal(i) = mean_signal(i) + mean(M(i, locs(j)-3:locs(j)+3))/length(locs);
    end
end
figure();
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 36, mean_signal, 'o', 'filled');
colorbar;
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end
hold on
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled') % selected dipole
legend('electrodes', 'selected dipole')
%% part f
Display_Potential_3D(9.2, mean_signal)
%% part g
alpha = 0.1;
% MNE
Q_MNE = G'*(G*G' + alpha*eye(21))^(-1)*mean_signal;
% WMNE
W_WMNE = kron(diag(Omega), eye(3));
Q_WMNE = (W_WMNE'*W_WMNE)^(-1)*G'*(G*(W_WMNE'*W_WMNE)^(-1)*G' + alpha*eye(N))^(-1)*mean_signal;
% LORETA
W_LORETA = kron(diag(Omega), eye(3)) * (B' * B) * kron(diag(Omega), eye(3));
Q_LORETA = (W_LORETA'*W_LORETA)^(-1)*G'*(G*(W_LORETA'*W_LORETA)^(-1)*G' + alpha*eye(N))^(-1)*mean_signal;
% SLORETA
Q_SLORETA = zeros(3*P, 1);
for i = 1:P
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end
%% part h
Q_MNE_max = -inf;
Q_MNE_max_idx = 0;
Q_WMNE_max = -inf;
Q_WMNE_max_idx = 0;
Q_LORETA_max = -inf;
Q_LORETA_max_idx = 0;
Q_SLORETA_max = -inf;
Q_SLORETA_max_idx = 0;
for i = 1:P
    if norm(Q_MNE(3*i-2:3*i)) > Q_MNE_max
        Q_MNE_max = norm(Q_MNE(3*i-2:3*i));
        Q_MNE_max_idx = i;
    end
    if norm(Q_WMNE(3*i-2:3*i)) > Q_WMNE_max
        Q_WMNE_max = norm(Q_WMNE(3*i-2:3*i));
        Q_WMNE_max_idx = i;
    end
    if norm(Q_LORETA(3*i-2:3*i)) > Q_LORETA_max
        Q_LORETA_max = norm(Q_LORETA(3*i-2:3*i));
        Q_LORETA_max_idx = i;
    end
    if norm(Q_SLORETA(3*i-2:3*i)) > Q_SLORETA_max
        Q_SLORETA_max = norm(Q_SLORETA(3*i-2:3*i));
        Q_SLORETA_max_idx = i;
    end
end
figure();
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); axis('equal');
xlabel('X', 'Interpreter','latex')
ylabel('Y', 'Interpreter','latex')
zlabel('Z', 'Interpreter','latex')
title('dipoles/electrodes positions', 'Interpreter','latex');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end

scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled')
scatter3(LocMat(1, Q_MNE_max_idx), LocMat(2, Q_MNE_max_idx), LocMat(3, Q_MNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_WMNE_max_idx), LocMat(2, Q_WMNE_max_idx), LocMat(3, Q_WMNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_LORETA_max_idx), LocMat(2, Q_LORETA_max_idx), LocMat(3, Q_LORETA_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_SLORETA_max_idx), LocMat(2, Q_SLORETA_max_idx), LocMat(3, Q_SLORETA_max_idx), 's', 'filled')
legend('dipoles', 'electrodes', 'selected dipole', 'predict (MNE)', 'predict (WMNE)', 'predict (LORETA)', 'predict (SLORETA)')

% location error
LocErr_MNE = mse(LocMat(:, n), LocMat(:, Q_MNE_max_idx));  
LocErr_WMNE = mse(LocMat(:, n), LocMat(:, Q_WMNE_max_idx));  
LocErr_LORETA = mse(LocMat(:, n), LocMat(:, Q_LORETA_max_idx)); 
LocErr_SLORETA = mse(LocMat(:, n), LocMat(:, Q_SLORETA_max_idx)); 
% angle error
e_q_MNE_predict = Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx) / norm(Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx));
AngErr_MNE = acos((e_q'*e_q_MNE_predict) / (norm(e_q)*norm(e_q_MNE_predict)))*180/pi;  
e_q_WMNE_predict = Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx) / norm(Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx));
AngErr_WMNE = acos((e_q'*e_q_WMNE_predict) / (norm(e_q)*norm(e_q_WMNE_predict)))*180/pi;  
e_q_LORETA_predict = Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx) / norm(Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx));
AngErr_LORETA = acos((e_q'*e_q_LORETA_predict) / (norm(e_q)*norm(e_q_LORETA_predict)))*180/pi;  
e_q_SLORETA_predict = Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx) / norm(Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx));
AngErr_SLORETA = acos((e_q'*e_q_SLORETA_predict) / (norm(e_q)*norm(e_q_SLORETA_predict)))*180/pi; 
%% part l Genetic Algorithm
x1 = [1289.0  -1.3882337679723866  -4.164700902593627  9.995282300223865];
q1 = x1(2:4);
e_q1 = q1 / norm(q1);
AngErr_GA = acos((eq_superficial'*e_q1') / (norm(e_q1')*norm(eq_superficial)))*180/pi; 
disp("AngErr_GA superficial dipole = "+AngErr_GA)

x2 = [853.0  6.23717949607848E-8  -5.373866997480836E-8  10.916849276991831];
q2 = x2(2:4);
e_q2 = q2/ norm(q2);
AngErr_GA_deep = acos((e_q_deep'*e_q2') / (norm(e_q2')*norm(e_q_deep)))*180/pi; 
disp("AngErr_GA deep dipole = "+AngErr_GA_deep)
%% part m
selected_dipoles = [1200:1204 1208:1212 1217:1221 1226:1230];

figure();
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); axis('equal');
xlabel('X', 'Interpreter','latex')
ylabel('Y', 'Interpreter','latex')
zlabel('Z', 'Interpreter','latex')
title('dipoles/electrodes positions', 'Interpreter','latex');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end
scatter3(LocMat(1, selected_dipoles), LocMat(2, selected_dipoles), LocMat(3, selected_dipoles),'MarkerFaceColor','yellow','MarkerEdgeColor','yellow') % selected dipole
e_q = zeros(3,length(selected_dipoles));

for i = 1:length(selected_dipoles)
    n = selected_dipoles(i);
    e_q(:,i) = LocMat(:,n)./norm(LocMat(:,n),2);
    n_eq = e_q(:,i);
    line([LocMat(1,n),LocMat(1,n)+n_eq(1)],[LocMat(2,n),LocMat(2,n)+n_eq(2)],[LocMat(3,n),LocMat(3,n)+n_eq(3)],'Color','yellow','linewidth',1.5);
end
legend('dipoles', 'electrodes', 'selected dipole')

M = zeros(21,10240);
for i = 1:length(selected_dipoles)
    n = selected_dipoles(i);
    q = Interictal(i,:);
    Q = e_q(:,i)*q;
    M = M + G(:,3*n-2:3*n)*Q;
end
disp_eeg(M, [], 256, ElecNames);
title('selected dipoles','interpreter','latex');

% find the channel with maximum voltage
max_amp = -inf;
max_idx = 0;
for i = 1:21
    if max(M(i, :) > max_amp)
        max_amp = max(M(i, :));
        max_idx = i;
    end
end
% locate windows with respect to channel with maximum peak voltage
[~, locs] = findpeaks(M(max_idx, :), 'MinPeakDistance', 10, 'MinPeakHeight', 70);
mean_signal = zeros(21,1);
for i = 1:21
    for j = 1:length(locs)
        mean_signal(i) = mean_signal(i) + mean(M(i, locs(j)-3:locs(j)+3))/length(locs);
    end
end
figure();
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 36, mean_signal, 'o', 'filled');
colorbar;
for i = 1:21
    text(ElecXYZ(1, i), ElecXYZ(2, i), ElecXYZ(3, i), ElecNames{i})
end
hold on
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled') % selected dipole
legend('electrodes', 'selected dipole')
%%
Display_Potential_3D(9.2, mean_signal)
%% part n
alpha = 0.1;
% MNE
Q_MNE = G'*(G*G' + alpha*eye(21))^(-1)*mean_signal;
% WMNE
Q_WMNE = (W_WMNE'*W_WMNE)^(-1)*G'*(G*(W_WMNE'*W_WMNE)^(-1)*G' + alpha*eye(N))^(-1)*mean_signal;
% LORETA
W_LORETA = kron(diag(Omega), eye(3)) * (B' * B) * kron(diag(Omega), eye(3));
Q_LORETA = (W_LORETA'*W_LORETA)^(-1)*G'*(G*(W_LORETA'*W_LORETA)^(-1)*G' + alpha*eye(N))^(-1)*mean_signal;
% SLORETA
Q_SLORETA = zeros(3*P, 1);
for i = 1:P
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

Q_MNE1 = zeros(1317, 1);
Q_WMNE1 = zeros(1317, 1);
Q_LORETA1 = zeros(1317, 1);
Q_SLORETA1 = zeros(1317, 1);
for i = 1:P
    Q_MNE1(i) = norm(Q_MNE(3*i-2:3*i));
    Q_WMNE1(i) = norm(Q_WMNE(3*i-2:3*i));
    Q_LORETA1(i) = norm(Q_LORETA(3*i-2:3*i));
    Q_SLORETA1(i) = norm(Q_SLORETA(3*i-2:3*i));
end
%% part o
trueDipoles = zeros(1317, 1);
trueDipoles(selected_dipoles) = 1;
[X,Y] = perfcurve(trueDipoles,Q_MNE1,1); 
figure();
plot(X, Y); xlabel('fpr'); ylabel('tpr');
hold on
%
[X,Y] = perfcurve(trueDipoles,Q_WMNE1,1); 
plot(X, Y);
%
[X,Y] = perfcurve(trueDipoles,Q_LORETA1,1); 
plot(X, Y);
%
[X,Y] = perfcurve(trueDipoles,Q_SLORETA1,1); 
plot(X, Y); 
legend('MNE', 'WMNE', 'LORETA', 'sLORETA', 'interpreter', 'latex')




















