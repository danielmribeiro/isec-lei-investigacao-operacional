%--------------------------------------------------------------------------
%                  Implementaçao do Metodo Simplex
%--------------------------------------------------------------------------
% Variaveis principais:
%--------------------------------------------------------------------------
% n = nº de variaveis originais
% m = nº de restriçoes funcionais
% A = matriz dos coeficientes técnicos
% b = vector dos termos independentes das restrições
% c = vector dos coeficientes das variaveis na FO
% x = vector com os indices de todas as variaveis
% xB = vector com os indices das variaveis básicas
% cB = vector com os coeficientes das variaveis básicas na FO
% Zjcj = vector com os valores da linha Zj-cj
% Z = valor da FO
% SBA = vector com os valores da SBA em cada iteração
%--------------------------------------------------------------------------

%limpa janela de comandos
clc
disp('-----------------------------------------------------------')
disp('        Resolução de um problema pelo método Simplex       ')
disp('-----------------------------------------------------------')
disp(' Assume-se que:                                            ')
disp(' -> Função objectivo está na forma de maximização          ')
disp(' -> Todas as restrições são de "<="                        ')
disp(' -> Todas as variáveis são >=0                             ')
disp('-----------------------------------------------------------')
% Le dados do problema
[n,m,c,A,b]=Le_dados;

%limpa novamente janela de comandos
clc

I=eye(m);       % Matriz identidade (mxm)
A=[A I];        % Matriz A do problema na forma aumentada
cs=zeros(1,m);  % Coeficientes das variaveis slack na FO
c=[c cs];       % Coeficientes de todas as variaveis na FO

xo=1:n;         % Indices das variaveis originais
xs=n+1:n+m;     % Indices das variaveis slack
x=[xo xs];      % Indices de todas as variaveis

xB=xs';         % xB e um vector coluna com os indices das VBs
cB=cs';         % cB e um vector coluna com os coeficientes das VBs na FO

SBA=[zeros(n,1);b];     % Inicializaçao da 1ª SBA (vector coluna)

Zjcj=zeros(1,n+m);      % Inicialização do vector Zjcj a zeros

termina=0;              % Controla a execuçao do ciclo
iteracao=1;             % Contabiliza o nº de iteraçoes
while ~termina
    
    % Completar código
    
    for j=1:m+n
        
        Zjcj(j)=cB'*A(:,j)-c(j);
    end
    
    Z=cB'*b;
    
    [valor_min coluna_pivot]=min(Zjcj);
    
    
    if valor_min >=0
        
        Apresenta_quadro_Simplex(n,m,c,xB,cB,A,b,Zjcj,Z,iteracao,0,0,0)
        termina=1;
        Apresenta_resultados_finais(n,m,SBA,Z)
        
        
    else
        
        for i=1:m
            if A(i,coluna_pivot)<=0
                q(i)=1E50;
                %q(i)=nan;
            else
                q(i)=b(i)/A(i,coluna_pivot);
            end
            
        end
        
        [val_min, linha_pivot]=min(q);
        
        
        Apresenta_quadro_Simplex(n,m,c,xB,cB,A,b,Zjcj,Z,iteracao,1,coluna_pivot,xB(linha_pivot))
        
        pause
        
        cB(linha_pivot)=c(coluna_pivot);
        xB(linha_pivot)=x(coluna_pivot);
        
        k=A(linha_pivot, coluna_pivot);
        
        A(linha_pivot,:)=A(linha_pivot,:)/k;
        b(linha_pivot)=b(linha_pivot)/k;
        
        
        for i=1:m
            
            if i~=linha_pivot
                k=A(i,coluna_pivot);
                A(i,:)=A(i,:)-k*A(linha_pivot,:);
                b(i)=b(i)-k*b(linha_pivot);
            end
        end
        
        SBA=zeros(m+n,1);
        
        for i=1:m
            SBA(xB(i))=b(i);
        end
        
        iteracao=iteracao+1;
        
                   
        
    end
    
  
    
end

