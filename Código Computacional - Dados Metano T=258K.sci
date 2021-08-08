//==============================================================================
//==============================================================================
//    Desenvolvimento de Código Computacional para Simulação de Formação   
//        de Hidratos de Gás com Base em Modelo Difusivo de Vlasov         
//==============================================================================
//==============================================================================
// Arquivo: Difusao_Hidrato_Ciclop.sci   
// Autor: Ivan Ramos
// Versão: 1.3
// Descrição: Difusão do Ciclopentano
// Status: Completo
//==============================================================================
//==============================================================================

clear;
clc;
format(25)

//==============================================================================
//______________________ 1. Definindo parâmetros _______________________________
//==============================================================================


//________________ 1.1 Parâmetros de Entrada do Usuário_________________________


//Condições Experimentais
R0 = 0.000026 ;        // Raio Inicial da Esfera de Gelo (m)
p = 6*(10^6) ;         // Pressão Experimental(Pa)
peq = 1.57*(10^6) ;         // Pressão (Pa)
T =  258 ;        //  Temperatura Experimental (K)
Teq = 258 ;         // Temperatura de Equilíbrio (K)
Z = 0.823088 ;         // fator de compressibilidade do gás hóspede nas condições experimentais
Zeq = 0.954117  ;       // fator de compressibilidade do gás hóspede no equilíbrio
R = 8.3144621 ;         //  Constante dos Gases (J/mol-K)
ttotal = 18000 ;           //  Tempo total de Experimento (min)

//Parâmetros dos Reagentes e do Hidrato
X = 7.8578*(10^3) ;         // densidade molar do hidrato de gás (mol/m^3)
W = 51*(10^3) ;         // densidade molar do gelo (mol/m^3)
Mw = 18.015 ;        // massa molar da água (g/mol)
Mh = 119.62625 ;        // massa molar do hidrato de gás (g/mol)
Mg =  16.04 ;        // Massa molar do gás (g/mol)
n =  5.75 ;        //  n (mols)
eh = 0.15 ;      // Porosidade do hidrato de gás

//Parâmetros da Formação
Deff = [9.4*(10^-17)] ;      // Coeficiente efetivo de difusão de gás no hidrato de gás (m^2/s)
K = [7*(10^-28)] ;         // taxa de formação de hidratos de gás a partir do gelo (m/mol-s)

//Pontos Experimentais
texpmin=60*[7.96296, 23.61111, 35.64815, 47.5, 61.85185, 79.81481, 98.98148, 116.11111, 137.22222, 161.48148, 182.77778, 211.75926] // Tempo dos pontos experimentais (min)
ntexp=[0.12222, 0.2254, 0.2746, 0.33206, 0.37683, 0.42571, 0.47556, 0.51651, 0.54952, 0.57778, 0.5981, 0.63206] // Valores de n(t) experimentais

//Parâmetros de Newton-Raphson
x0 = 0.00001;
tol = 10^(-10);


//_______________ 1.2 Cálculo dos demais parâmetros iniciais____________________


//Definição da escala do gráfico para melhor visualização e cálculo dos resultados
thms=2
if ttotal >= 300 then
    ttotal = ttotal/60
    thms=1
    disp("tempo em minutos muito grande, gráfico será plotado em horas")
elseif ttotal <= 20 then
    ttotal = ttotal*60
    thms=3
    disp('tempo em minutos muito curto, gráfico será plotado em segundos')
end

//Parâmetros Independentes do tempo
deltaform = ((p/(Z*R*T))- (peq/(Zeq*R*Teq))) ;  //quantidade volumétrica de hidrato formado corrigida para o caso do experimento com ciclopentano
psi = ((W*Mw)/(X*Mh*(1-eh)))*(1 + (Mg/(n*Mw))) ;
A = (Deff)/(K*(W^n)) ;
B = 1/(2*(psi-1)) ;
C = (R0^3)*psi ;


//==============================================================================
//___________________________2. Criando funções ________________________________
//==============================================================================

//________________________________2.1[Função f(d)]______________________________
//[Função do numerador do algoritmo de Newton-Raphson]

function y=f(x)
y = 2*B*(C - (x^3)*(psi-1))^(2/3) + x^2 + 2*A*x - E;
endfunction

//_____________________________2.2[Função F'(d)]________________________________
//[Função do denominador do algoritmo de Newton-Raphson]

function y=df(x)
y = 2*(A + x - (((2*B*(x^2))*(psi-1))/((C - (x^3)*(psi-1))^(1/3))));
endfunction

//___________________2.3 [Função Método Newton-Raphson]_________________________
//[]

function ksi=newtonraphson(x0, tol)
i=0;
ksi = x0 - (f(x0)/df(x0));

while (abs(ksi-x0)>tol)
x0=ksi;
ksi = x0 - (f(x0)/df(x0));
i=i+1;
end

endfunction

function y=raiodegelo(t)
    F = (Deff*deltaform*t)/(X*(1-eh)) ;
    E = (R0*(2*A + R0*(2*B +1)) - 2*F) ;
    y=newtonraphson(x0, tol);
endfunction

function y=graudeformacao(t)
    y = 1 - (raiodegelo(t)/R0)^3
endfunction

function y=raioexternodacamada(t)
    y = ((psi*((R0^3)-(raiodegelo(t)^3)))+(raiodegelo(t)^3))^(1/3)
endfunction

//==============================================================================
//_______________________________3. Cálculo ____________________________________
//==============================================================================

tk=0
i=1
while tk <= ttotal
    
    if thms==1
        tseg= tk*60*60 ;
    elseif thms==3
        tseg= tk
    else
        tseg= tk*60
    end 
    
    t(i)=tk

    ksi(i,1)= tk
    ksi(i,2)=raiodegelo(tseg)  //armazenando valores do raio de gelo em um vetor
    
    nt(i,1)= tk
    nt(i,2)= graudeformacao(tseg) // calculo do grau de formação
    
    rt(i,1)= tk
    rt(i,2)= raioexternodacamada(tseg) //((psi*((R0^3)-(ksi(i,2)^3)))+(ksi(i,2)^3))^(1/3) //calculo do raio da esfera
    
    i=i+1
    if thms==1
        tk=tk+0.1
    else
        tk=tk+1
    end
end



//==============================================================================
//______________________3. Ferramentas de Gráfico ______________________________
//==============================================================================
//[Ferramentas para plotagem dos gráficos]


if thms==1 then
    texpmin=texpmin/60
end
if thms==3 then
    texpmin=texpmin*60
end

i=1:size(t, 1)

//ferramenta gráfica
plot(t(i), nt(i,2), 'k-', 'LineWidth', 2);

//plotagem dos valores experimentais
scatter(texpmin, ntexp, msizes = 20, "red")

a=gca();
a.data_bounds = [0,0;ttotal,1]
