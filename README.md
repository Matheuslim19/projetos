# projetos
%parametros
Lx = 10;                % Comprimento em x
Ly = 2;                 % Comprimento em y
t_final = 48*30;           % Tempo final 1 Semana (7 dias)
dt = 0.01;              %passo do tempo
C0 = 10;                % Concentração inicial
u = 0.2;                % Velocidade em x 
v = 0.01;               % Velocidade em y
kx = 0.00002;           % Coeficiente de difusão em x
ky = 0.00001;           % Coeficiente de difusão em y


% Parâmetros do método das diferenças finitas
Nx = 100;               % Número de pontos em x
Ny = 20;                % Número de pontos em y
dx = Lx/(Nx-1);         % Intervalo de espaço em x
dy = Ly/(Ny-1);         % Intervalo de espaço em y


% método das diferenças finitas
alpha_x = kx*dt/(dx^2); % Parâmetro de difusão em x
alpha_y = ky*dt/(dy^2); % Parâmetro de difusão em y
beta_x = ((u*dt)/dx);
beta_y = ((v*dt)/dy);


% CONDIÇÃO DE CONTORNO
C = zeros(Nx,Ny);       % Matriz de concentrações
C(:,1) = 0;             % Condição de contorno em y=0
C(:,end) = 0;           % Condição de contorno em y=Ly
C(1,:) = 0;             % Condição de contorno em x=0
C(end,:) = 0;           % Condição de contorno em x=Lx
C(1, (Ny/2)) = C0;          % Condição inicial com ponto de derramamento

t = 0;
while t < t_final
    
    % Atualização do tempo
    t = t + dt;
    
% Cálculo das derivadas parciais
% Derivada em x
dCdx = diff(C, 1, 1) / dx;
dCdx = padarray(dCdx, [1, 0]); % preencher a última linha com zeros

% Derivada em y
dCdy = diff(C, 1, 2) / dy;
dCdy = [dCdy, zeros(size(dCdy, 1), 1)]; % preencher a última coluna com zeros
    
    
 % Cálculo do termo difusivo com diferença central
    d2Cdx2 = zeros(Nx,Ny);
    d2Cdy2 = zeros(Nx,Ny);
    
    for i = 2:Nx-1
        d2Cdx2(i,:) = (C(i+1,:)-2*C(i,:)+C(i-1,:))/dx^2;
    end
    d2Cdx2(:,1) = 0;
    d2Cdx2(:,end) = 0;

    for j = 2:Ny-1
        d2Cdy2(:,j) = (C(:,j+1)-2*C(:,j)+C(:,j-1))/dy^2;
    end
    d2Cdy2(1,:) = 0;
    d2Cdy2(end,:) = 0;

  
  
    % Cálculo do termo advectivo Método UP-Wind de 1° ordem
    adv_x = zeros(Nx,Ny);
    adv_y = zeros(Nx,Ny);
    
    for i = 2:Nx-1
        if u>=0
            adv_x(i,:) = -u*(C(i,:)-C(i-1,:))/dx;
        else
            adv_x(i,:) = -u*(C(i+1,:)-C(i,:))/dx;
        end
    end

    for j = 2:Ny-1
        if v>=0
            adv_y(:,j) = -v*(C(:,j)-C(:,j-1))/dy;
        else
            adv_y(:,j) = -v*(C(:,j+1)-C(:,j))/dy;
        end
    end
   
    % atualização da concentração( Método de Euler explícito )
    C = C + alpha_x*(d2Cdx2 + adv_x) + alpha_y*(d2Cdy2 + adv_y) * dt;
    
end

% Gerando a malha de plotagem
[x, y] = meshgrid(0:dx:Lx, 0:dy:Ly);
% Plot da concentração
figure(1)
contourf(x, y, C', 50, 'LineStyle', 'none')
colormap jet
colorbar
xlabel('x (Km)')
ylabel('y (Km)')
title('Concentração de óleo no tempo')
hold on
grid on
% Definindo o nível do contorno
contour(x, y, C', [0.8 1 1.2], 'LineColor', 'w', 'LineWidth', 1.5)
