clc
clear
close all

N = input("Enter N: ");
if N > 1
    e1 = input("Enter e1: ");
    e2 = input("Enter e2: ");
    dc = input("Do you want the design to be DC-drop-free? Press 1 for yes, 0 for no.");
else
    fprintf("N must be equal or greater than 2")
end

e12 = e1*e2;

syms w
syms s

even = 0;
if rem(N,2) == 0
    even = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% whit DC-drop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (dc == 0 || even == 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F11 & F12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a11 = sinh(1/N*asinh(1/e1));
    b11 = cosh(1/N*asinh(1/e1));
    
    a12 = sinh(1/N*asinh(1/e2));
    b12 = cosh(1/N*asinh(1/e2));

    pk11 = zeros(1, N);
    pk12 = zeros(1, N);
    for i = 1:N
        pk11(i) = -a11*sin((2*i-1)*pi/(2*N)) + 1j*b11*cos((2*i-1)*pi/(2*N));
        pk12(i) = -a12*sin((2*i-1)*pi/(2*N)) + 1j*b12*cos((2*i-1)*pi/(2*N));
    end

    NUM11 = 1/(e1*2^(N-1));
    NUM12 = 1/(e2*2^(N-1));

    DEN11 = [1 -pk11(1)];
    DEN12 = [1 -pk12(1)];
    for i = 2:N
        DEN11_temp = [1 -pk11(i)];
        DEN12_temp = [1 -pk12(i)];
        DEN11 = conv(DEN11, DEN11_temp);
        DEN12 = conv(DEN12, DEN12_temp);
    end
    
    NUM11
    NUM12
    DEN11
    DEN12
    NUM11(imag(NUM11)<1e-10) = real(NUM11);
    NUM12(imag(NUM12)<1e-10) = real(NUM12);
    DEN11(imag(DEN11)<1e-10) = real(DEN11);
    DEN12(imag(DEN12)<1e-10) = real(DEN12);
    NUM11
    NUM12
    DEN11
    DEN12
    TF_F11 = tf(NUM11, DEN11);
    TF_F12 = tf(NUM12, DEN12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NUM1 = conv(NUM11, NUM12);
    DEN1 = conv(DEN11, DEN12);
    TF_F1 = tf(NUM1, DEN1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a22 = sinh(1/(2*N)*asinh(1/e12));
    b22 = cosh(1/(2*N)*asinh(1/e12));

    pk2 = zeros(1, 2*N);

    for i = 1:2*N 
        pk2(i) = -a22*sin((2*i-1)*pi/(2*2*N)) + 1j*b22*cos((2*i-1)*pi/(2*2*N));  
    end

    NUM2 = 1/(e12*2^(2*N-1));
    DEN2 = [1 -pk2(1)];
    for i = 2:2*N
        DEN2_temp = [1 -pk2(i)];
        DEN2 = conv(DEN2, DEN2_temp);
    end

    NUM2(imag(NUM2)<1e-10) = real(NUM2);
    DEN2(imag(DEN2)<1e-10) = real(DEN2);
    TF_F2 = tf(NUM2, DEN2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLES PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1)
    
    subplot(2,1,1)
    for i = 1:N
        title("poles of F1")
        xlabel("Real")
        ylabel("Imaginary")
        plot(pk11(i), 'x')    
        hold on
        grid on
        plot(pk12(i), 'x')
        hold on
    end
    
    subplot(2,1,2)
    for i = 1:2*N
        title("poles of F2")
        xlabel("Real")
        ylabel("Imaginary")
        plot(pk2(i), 'x')
        grid on
        hold on
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELAY PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(2)

    subplot(2,1,1)
    [h_F1, w_F1] = freqs(DEN1, NUM1, 1000);
    grpdel_F1 = diff(unwrap(angle(h_F1)))./diff(w_F1);
    semilogx(w_F1(2:end), grpdel_F1)
    xlabel('Frequency (rad/s)')
    ylabel('Group delay (s)')
    title('F1')
    grid on
    
    subplot(2,1,2)    
    [h_F2, w_F2] = freqs(DEN2, NUM2, 1000);
    grpdel_F2 = diff(unwrap(angle(h_F2)))./diff(w_F2);
    semilogx(w_F2(2:end), grpdel_F1)
    xlabel('Frequency (rad/s)')
    ylabel('Group delay (s)')
    title('F2')
    grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAGNITUDE CHARACTERISTIC PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(3)

    subplot(2,1,1)
    DEN1_w = poly2sym(DEN1, w);
    DEN1_jw = subs(DEN1_w, w, sqrt(-1)*w);
    F1_jw = NUM1/DEN1_jw;
    title("Magnitude Characteristic of F1")
    xlabel("w")
    ylabel("|F1(jw)|")
    grid on
    fplot(abs(F1_jw))

    subplot(2,1,2)
    DEN2_w = poly2sym(DEN2, w);
    DEN2_jw = subs(DEN2_w, w, sqrt(-1)*w);
    F2_jw = NUM2/DEN2_jw;
    title("Magnitude Characteristic of F2")
    xlabel("w")
    ylabel("|F2(jw)|")
    grid on
    fplot(abs(F2_jw))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% without DC-drop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif (dc == 1 && even == 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T(N)(w) function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TNw_N = sym('T1',[1 N+1]);
TNw_N(3) = 2*w^2 - 1;
TNw_N(2) = w;
TNw_N(1) = 1;

if (N == 3)
    TNw_N(4) = 4*w^3 - 3*w;
elseif (N > 3)
    TNw_N(4) = 4*w^3 - 3*w;

    TNw_1_2ndpower_arrays = floor(log(N)/log(2));
    for k = 2:TNw_1_2ndpower_arrays
        TNw_N(2^k + 1) = 2*(TNw_N(2*(k-1) + 1))^2 - 1;
    end
if (N > 4)
    for i = 5:(N+1)
        if (log(i-1)/log(2)) - floor(log(i-1)/log(2)) ~= 0
            TNw_N(i) = 2*w*TNw_N(i-1) - TNw_N(i-2);
        end
    end
end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T(2N)(w) function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TNw_2N = sym('T2', [1 2*N+1]);
    TNw_2N(4) = 4*w^3 - 3*w;
    TNw_2N(3) = 2*w^2 - 1;
    TNw_2N(2) = w;
    TNw_2N(1) = 1;
    
    TNw_2N_2ndpower_arrays = floor(log(2*N)/log(2));
        for i = 2:TNw_2N_2ndpower_arrays
            TNw_2N(2^i + 1) = 2*(TNw_2N(2*(i-1) + 1))^2 - 1;
        end
    if (N > 4)
        for i = 5:(2*N+1)
            if (log(2*i-1)/log(2)) - floor(log(2*i-1)/log(2)) ~= 0
                TNw_2N(i) = 2*w*TNw_2N(i-1) - TNw_2N(i-2);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% |F11|^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    LNw_N = subs(TNw_N(N+1), w, sqrt( (cos(pi/(2*N))*w)^2 + (sin(pi/(2*N)))^2 ));
    F11_squared_wo_DC_drop_DEN = 1 + (e1*LNw_N)^2;
    F11_squared_wo_DC_drop_DEN_s = subs(F11_squared_wo_DC_drop_DEN, w, s/sqrt(-1));
    poles_F11_squared_wo_DC_drop = double(solve(F11_squared_wo_DC_drop_DEN_s == 0 , s));
    poles_F11_squared_wo_DC_drop = unique(poles_F11_squared_wo_DC_drop);

    poles_F11_wo_DC_drop = zeros(1, N);
    K11_wo_DC_drop = 1;     
    j = 1;
    for i = 1:2*N
        if real(poles_F11_squared_wo_DC_drop(i)) < 0
            poles_F11_wo_DC_drop(j) = poles_F11_squared_wo_DC_drop(i);
            K11_wo_DC_drop = K11_wo_DC_drop*(0-poles_F11_wo_DC_drop(j));
            j = j+1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% |F12|^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F12_squared_wo_DC_drop_DEN = 1 + (e2*LNw_N)^2;
    F12_squared_wo_DC_drop_DEN_s = subs(F12_squared_wo_DC_drop_DEN, w, s/sqrt(-1));
    poles_F12_squared_wo_DC_drop = double(solve(F12_squared_wo_DC_drop_DEN_s == 0 , s));
    poles_F12_squared_wo_DC_drop = unique(poles_F12_squared_wo_DC_drop);

    poles_F12_wo_DC_drop = zeros(1, N);
    K12_wo_DC_drop = 1;     
    j = 1;
    for i = 1:2*N
        if real(poles_F12_squared_wo_DC_drop(i)) < 0
            poles_F12_wo_DC_drop(j) = poles_F12_squared_wo_DC_drop(i);
            K12_wo_DC_drop = K12_wo_DC_drop*(0-poles_F12_wo_DC_drop(j));
            j = j + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F11 and F12 without DC drop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NUM11_wo_DC_drop = K11_wo_DC_drop;
    NUM12_wo_DC_drop = K12_wo_DC_drop;

    DEN11_wo_DC_drop = [1 -poles_F11_wo_DC_drop(1)];
    DEN12_wo_DC_drop = [1 -poles_F12_wo_DC_drop(1)];
    for i = 2:N
        DEN11_wo_DC_drop_temp = [1 -poles_F11_wo_DC_drop(i)];
        DEN12_wo_DC_drop_temp = [1 -poles_F12_wo_DC_drop(i)];
        DEN11_wo_DC_drop = conv(DEN11_wo_DC_drop, DEN11_wo_DC_drop_temp);
        DEN12_wo_DC_drop = conv(DEN12_wo_DC_drop, DEN12_wo_DC_drop_temp);
    end

    NUM11_wo_DC_drop(imag(NUM11_wo_DC_drop)<1e-10) = real(NUM11_wo_DC_drop);
    NUM12_wo_DC_drop(imag(NUM12_wo_DC_drop)<1e-10) = real(NUM12_wo_DC_drop);
    DEN11_wo_DC_drop(imag(DEN11_wo_DC_drop)<1e-10) = real(DEN11_wo_DC_drop);
    DEN12_wo_DC_drop(imag(DEN12_wo_DC_drop)<1e-10) = real(DEN12_wo_DC_drop);

    TF_F11_wo_DC_drop = tf(NUM11_wo_DC_drop, DEN11_wo_DC_drop);
    TF_F12_wo_DC_drop = tf(NUM12_wo_DC_drop, DEN12_wo_DC_drop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NUM1 = conv(NUM11_wo_DC_drop, NUM12_wo_DC_drop);
    DEN1 = conv(DEN11_wo_DC_drop, DEN12_wo_DC_drop);
    TF_F1_wo_DC_drop = tf(NUM1, DEN1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% |F2|^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    LNw_2N = subs(TNw_2N(2*N+1), w, sqrt( (cos(pi/(2*2*N))*w)^2 + (sin(pi/(2*2*N)))^2 ));
    F2_squared_wo_DC_drop_DEN = 1 + (e12*LNw_2N)^2;
    F2_squared_wo_DC_drop_DEN_s = subs(F2_squared_wo_DC_drop_DEN, w, s/sqrt(-1));
    poles_F2_squared_wo_DC_drop = double(solve(F2_squared_wo_DC_drop_DEN_s == 0 , s));
    poles_F2_squared_wo_DC_drop = unique(poles_F2_squared_wo_DC_drop);

    poles_F2_wo_DC_drop = ones(1, 2*N); 
    K2_wo_DC_drop = 1;
    j = 1;
    for k = 1:2*2*N
        if real(poles_F2_squared_wo_DC_drop(k)) < 0
            poles_F2_wo_DC_drop(j) = poles_F2_squared_wo_DC_drop(k);
            K2_wo_DC_drop = K2_wo_DC_drop*(0-poles_F2_wo_DC_drop(j));
            j = j + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F2 without DC drop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NUM2 =  K2_wo_DC_drop;

    DEN2 = [1 -poles_F2_wo_DC_drop(1)];
    for i = 2:2*N
        DEN2_wo_DC_drop_temp = [1 -poles_F2_wo_DC_drop(i)];
        DEN2 = conv(DEN2, DEN2_wo_DC_drop_temp);
    end

    NUM2(imag(NUM2)<1e-10) = real(NUM2);
    DEN2(imag(DEN2)<1e-10) = real(DEN2);

    TF_F2_wo_DC_drop = tf(NUM2, DEN2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLES PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1)

    subplot(2,1,1)
    for i = 1:N
        title("poles of F1 without DC drop")
        xlabel("Real")
        ylabel("Imaginary")
        plot(poles_F11_wo_DC_drop(i), 'x')    
        hold on
        grid on
        plot(poles_F12_wo_DC_drop(i), 'x')
        hold on
    end
    
    subplot(2,1,2)
    for i = 1:2*N
        title("poles of F2 without DC drop")
        xlabel("Real")
        ylabel("Imaginary")
        grid on
        plot(poles_F2_wo_DC_drop(i), 'x')
        hold on
    end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELAY PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(2)

    subplot(2,1,1)
    [h_F1_wo_drop, w_F1_wo_drop] = freqs(DEN1, NUM1, 1000);
    grpdel_F1_wo_drop = diff(unwrap(angle(h_F1_wo_drop)))./diff(w_F1_wo_drop);
    semilogx(w_F1_wo_drop(2:end), grpdel_F1_wo_drop)
    xlabel('Frequency (rad/s)')
    ylabel('Group delay (s)')
    grid on
    title('F1')
    
    subplot(2,1,2)    
    [h_F2_wo_drop, w_F2_wo_drop] = freqs(DEN2, NUM2, 1000);
    grpdel_F2_wo_drop = diff(unwrap(angle(h_F2_wo_drop)))./diff(w_F2_wo_drop);
    semilogx(w_F2_wo_drop(2:end), grpdel_F2_wo_drop)
    xlabel('Frequency (rad/s)')
    ylabel('Group delay (s)')
    grid on
    title('F2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAGNITUDE CHARACTERISTIC PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(3)

    subplot(2,1,1)
    DEN1_w_wo_DC_drop = poly2sym(DEN1, w);
    DEN1_jw_wo_DC_drop = subs(DEN1_w_wo_DC_drop, w, sqrt(-1)*w);
    F1_jw_wo_DC_drop = NUM1/DEN1_jw_wo_DC_drop;
    title("Magnitude Characteristic of F1 without DC drop")
    xlabel("w")
    ylabel("|F1(jw)|")
    grid on
    fplot(abs(F1_jw_wo_DC_drop))

    subplot(2,1,2)
    DEN2_w_wo_DC_drop = poly2sym(DEN2, w);
    DEN2_jw_wo_DC_drop = subs(DEN2_w_wo_DC_drop, w, sqrt(-1)*w);
    F2_jw_wo_DC_drop = NUM2/DEN2_jw_wo_DC_drop;
    title("Magnitude Characteristic of F2 without DC drop")
    xlabel("w")
    ylabel("|F2(jw)|")
    grid on
    fplot(abs(F2_jw_wo_DC_drop))    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rg_denorm = input("Enter Rg: ");
RL_denorm = input("Enter RL: ");
Rg = Rg_denorm/RL_denorm;
RL = 1;

t1_s = 2*sqrt(Rg)/(1+Rg)*NUM1/poly2sym(DEN1, s);

p1_s_squared = 1 - t1_s*subs(t1_s, s, -s);
p1_s_zeros = double(solve(p1_s_squared == 0));
p1_s_zeros1 = p1_s_zeros(real(p1_s_zeros)<=0);

p1_s_zeros2 = zeros(1, 2*N);
j = 1;
for i=1:size(p1_s_zeros1,1)-2
    if real(p1_s_zeros1(i))==0
        if abs(imag(p1_s_zeros1(i))-imag(p1_s_zeros1(i+2)))<1e-5
            p1_s_zeros2(j) = p1_s_zeros1(i);
            j = j + 1;
        end
    end
end
for i=1:size(p1_s_zeros1,1)
    if real(p1_s_zeros1(i))~=0
        p1_s_zeros2(j) = p1_s_zeros1(i);
        j = j + 1;
    end
end

p1_0 = (1-Rg)/(1+Rg);
if (p1_0 == 0)
    fprintf("circuit of F1 has dual \r");
    K_p1 = sqrt(limit(p1_s_squared, s, inf));
elseif (p1_0 > 0)
    K_p1 = sqrt(limit(p1_s_squared, s, inf));
elseif (p1_0 < 0)
    K_p1 = -sqrt(limit(p1_s_squared, s, inf));
end

p1_NUM = [1 -p1_s_zeros2(1)];
[p1_s_zeros2_rows, p1_s_zeros2_columns] = size(p1_s_zeros2);
for i = 2:p1_s_zeros2_columns
    p1_NUM_temp = [1 -p1_s_zeros2(i)];
    p1_NUM = conv(p1_NUM, p1_NUM_temp);
end
p1_NUM(imag(p1_NUM)<1e-10) = real(p1_NUM);

p1_NUM_s = poly2sym(p1_NUM, s);
p1 = K_p1*p1_NUM_s/poly2sym(DEN1, s);

A1_B1 = sqrt(Rg)*(1+p1)/t1_s;
C1_D1 = sqrt(Rg)*(1-p1)/t1_s;

A1_B1_odd_even = sym2poly(simplify(A1_B1));
C1_D1_odd_even = sym2poly(simplify(C1_D1));

[A1_B1_odd_even_row, A1_B1_odd_even_column] = size(A1_B1_odd_even);
[C1_D1_odd_even_row, C1_D1_odd_even_column] = size(C1_D1_odd_even);

if rem(A1_B1_odd_even_column,2) == 1
    j = 1;
    A1_B1_even = zeros(1, N+1);
    for i=1:2:A1_B1_odd_even_column
        A1_B1_even(j) = A1_B1_odd_even(i);
        j = j + 1;
    end
    j = 1;
    A1_B1_odd = zeros(1, N);
    for i=2:2:A1_B1_odd_even_column
        A1_B1_odd(j) = A1_B1_odd_even(i);
        j = j + 1;
    end
elseif rem(A1_B1_odd_even_column,2) == 0
    j = 1;
    A1_B1_even = zeros(1, N);
    for i=1:2:A1_B1_odd_even_column
        A1_B1_even(j) = A1_B1_odd_even(i);
        j = j + 1;
    end
    j = 1;
    A1_B1_odd = zeros(1, N+1);
    for i=2:2:A1_B1_odd_even_column
        A1_B1_odd(j) = A1_B1_odd_even(i);
        j = j + 1;
    end
end

if rem(C1_D1_odd_even_column,2) == 1
    j = 1;
    C1_D1_even = zeros(1, N);
    for i=1:2:C1_D1_odd_even_column
        C1_D1_even(j) = C1_D1_odd_even(i);
        j = j + 1;
    end
    j = 1;
    C1_D1_odd = zeros(1, N+1);
    for i=2:2:C1_D1_odd_even_column
        C1_D1_odd(j) = C1_D1_odd_even(i);
        j = j + 1;
    end
elseif rem(C1_D1_odd_even_column,2) == 0
    j = 1;
    C1_D1_even = zeros(1, N);
    for i=2:2:C1_D1_odd_even_column
        C1_D1_even(j) = C1_D1_odd_even(i);
        j = j + 1;
    end
    j = 1;
    C1_D1_odd = zeros(1, N);
    for i=1:2:C1_D1_odd_even_column
        C1_D1_odd(j) = C1_D1_odd_even(i);
        j = j + 1;
    end
end

A1_temp = poly2sym(A1_B1_even, s);
A1 = subs(A1_temp, s, s^2);
B1_temp = poly2sym(A1_B1_odd, s);
B1 = s*subs(B1_temp, s, s^2);
C1_temp = poly2sym(C1_D1_odd, s);
C1 = s*subs(C1_temp, s, s^2);
D1_temp = poly2sym(C1_D1_even, s);
D1 = subs(D1_temp, s, s^2);

T1 = vpa([A1 B1; C1 D1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t2_s = 2*sqrt(Rg)/(1+Rg)*NUM2/poly2sym(DEN2, s);

p2_s_squared = 1 - t2_s*subs(t2_s, s, -s);
p2_s_zeros = double(solve(p2_s_squared == 0));
p2_s_zeros1 = p2_s_zeros(real(p2_s_zeros)<=1e-10);

p2_s_zeros2 = zeros(1, 2*N);
j = 1;
for i=1:size(p2_s_zeros1,1)-2
    if real(p2_s_zeros1(i))==0
        if abs(imag(p2_s_zeros1(i))-imag(p2_s_zeros1(i+2)))<1e-5
            p2_s_zeros2(j) = p2_s_zeros1(i);
            j = j + 1;
        end
    end
end
for i=1:size(p2_s_zeros1,1)
    if real(p2_s_zeros1(i))~=0
        p2_s_zeros2(j) = p2_s_zeros1(i);
        j = j + 1;
    end
end

p2_0 = (1-Rg)/(1+Rg);
if (p2_0 == 0)
    fprintf("circuit of F2 has dual \r");
    K_p2 = sqrt(limit(p2_s_squared, s, inf));
elseif (p2_0 > 0)
    K_p2 = sqrt(limit(p2_s_squared, s, inf));
elseif (p2_0 < 0)
    K_p2 = -sqrt(limit(p2_s_squared, s, inf));
end

p2_NUM = [1 -p2_s_zeros2(1)];
[p2_s_zeros2_rows, p2_s_zeros2_columns] = size(p2_s_zeros2);
for i = 2:p2_s_zeros2_columns
    p2_NUM_temp = [1 -p2_s_zeros2(i)];
    p2_NUM = conv(p2_NUM, p2_NUM_temp);
end
p2_NUM(imag(p2_NUM)<1e-10) = real(p2_NUM);

p2_NUM_s = poly2sym(p2_NUM, s);
p2 = K_p2*p2_NUM_s/poly2sym(DEN2, s);

A2_B2 = sqrt(Rg)*(1+p2)/t2_s;
C2_D2 = sqrt(Rg)*(1-p2)/t2_s;

A2_B2_odd_even = sym2poly(simplify(A2_B2));
C2_D2_odd_even = sym2poly(simplify(C2_D2));

[A2_B2_odd_even_row, A2_B2_odd_even_column] = size(A2_B2_odd_even);
[C2_D2_odd_even_row, C2_D2_odd_even_column] = size(C2_D2_odd_even);

if rem(A2_B2_odd_even_column,2) == 1
    j = 1;
    A2_B2_even = zeros(1, N+1);
    for i=1:2:A2_B2_odd_even_column
        A2_B2_even(j) = A2_B2_odd_even(i);
        j = j + 1;
    end
    j = 1;
    A2_B2_odd = zeros(1, N);
    for i=2:2:A2_B2_odd_even_column
        A2_B2_odd(j) = A2_B2_odd_even(i);
        j = j + 1;
    end
elseif rem(A2_B2_odd_even_column,2) == 0
    j = 1;
    A2_B2_odd = zeros(1, N);
    for i=1:2:A2_B2_odd_even_column
        A2_B2_odd(j) = A2_B2_odd_even(i);
        j = j + 1;
    end
    j = 1;
    A2_B2_even = zeros(1, N+1);
    for i=2:2:A2_B2_odd_even_column
        A2_B2_even(j) = A2_B2_odd_even(i);
        j = j + 1;
    end
end

if rem(C2_D2_odd_even_column,2) == 1
    j = 1;
    C2_D2_even = zeros(1, N+1);
    for i=1:2:C2_D2_odd_even_column
        C2_D2_even(j) = C2_D2_odd_even(i);
        j = j + 1;
    end
    j = 1;
    C2_D2_odd = zeros(1, N);
    for i=2:2:C2_D2_odd_even_column
        C2_D2_odd(j) = C2_D2_odd_even(i);
        j = j + 1;
    end
elseif rem(C2_D2_odd_even_column,2) == 0
    j = 1;
    C2_D2_even = zeros(1, N);
    for i=2:2:C2_D2_odd_even_column
        C2_D2_even(j) = C2_D2_odd_even(i);
        j = j + 1;
    end
    j = 1;
    C2_D2_odd = zeros(1, N);
    for i=1:2:C2_D2_odd_even_column
        C2_D2_odd(j) = C2_D2_odd_even(i);
        j = j + 1;
    end
end

A2_temp = poly2sym(A2_B2_even, s);
A2 = subs(A2_temp, s, s^2);
B2_temp = poly2sym(A2_B2_odd, s);
B2 = s*subs(B2_temp, s, s^2);
C2_temp = poly2sym(C2_D2_odd, s);
C2 = s*subs(C2_temp, s, s^2);
D2_temp = poly2sym(C2_D2_even, s);
D2 = subs(D2_temp, s, s^2);

T2 = vpa([A2 B2; C2 D2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% decomposition of T1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1_matrix = cell(1, 2*N);
b1_T = sym('b1', [1 2*N]);
b1_T_final = zeros(1, 2*N);
e1_T = sym('e1', [1 2*N]);
e1_T_final = zeros(1, 2*N);



T1_temp = T1;
A1_inv_T = [1 -b1_T(1)*s; -e1_T(1)*s 1]*T1_temp;
for i = 1:2*N
    A1_eq = cell(1, 4);
    A1_eq_row = zeros(1, 4);
    A1_eq_column = zeros(1, 4);
    for j = 1:4
        A1_eq{j} = coeffs(A1_inv_T(j), s, 'ALL');
        [A1_eq_row(j), A1_eq_column(j)] = size(A1_eq{j});
        if A1_eq_column(j) > 2*N-(i-1)
            if j == 1
                if i == 1
                    b1_T_final(i) = double(solve(A1_eq{1}(1) == 0));
                elseif i > 1
                    if e1_T_final(i-1) > 1e-7
                        for k = 1:2:A1_eq_column(j)
                            if double(solve(A1_eq{j}(k) == 0, b1_T(i))) > 1e-7
                                b1_T_final(i) = double(solve(A1_eq{j}(k) == 0, b1_T(i)));
                                break
                            end
                        end
                    end
                end
                
            elseif j == 2 || j == 4
                if i == 1
                    if j == 4
                        e1_T_final(i) = double(solve(A1_eq{4}(1) == 0));
                    end
                elseif i > 1
                    if b1_T_final(i-1) > 1e-7
                        for k = 1:2:A1_eq_column(j)
                            if double(solve(A1_eq{j}(k) == 0, e1_T(i))) > 1e-7
                                e1_T_final(i) = double(solve(A1_eq{j}(k) == 0, e1_T(i)));
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    A1_matrix{i} = [1 b1_T_final(i)*s; e1_T_final(i)*s 1];
    T1_temp = [1 -b1_T_final(i)*s; -e1_T_final(i)*s 1]*T1_temp;
    if i < 2*N
        A1_inv_T = [1 -b1_T(i+1)*s; -e1_T(i+1)*s 1]*T1_temp;
    end
end

b1_T_final
e1_T_final

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% decomposition of T2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A2_matrix = cell(1, 2*N);
b2_T = sym('b1', [1 2*N]);
b2_T_final = zeros(1, 2*N);
e2_T = sym('e1', [1 2*N]);
e2_T_final = zeros(1, 2*N);

T2_temp = T2;
A2_inv_T = [1 -b2_T(1)*s; -e2_T(1)*s 1]*T2_temp;
for i = 1:2*N
    A2_eq = cell(1, 4);
    A2_eq_row = zeros(1, 4);
    A2_eq_column = zeros(1, 4);
    for j = 1:4
        A2_eq{j} = coeffs(A2_inv_T(j), s, 'ALL');
        [A2_eq_row(j), A2_eq_column(j)] = size(A2_eq{j});
        if A2_eq_column(j) > 2*N-(i-1)
            if j == 1
                if i == 1
                    b2_T_final(i) = double(solve(A2_eq{1}(1) == 0));
                elseif i > 1
                    if e2_T_final(i-1) > 1e-7
                        for k = 1:2:A2_eq_column(j)
                            if double(solve(A2_eq{j}(k) == 0, b2_T(i))) > 1e-7
                                b2_T_final(i) = double(solve(A2_eq{j}(k) == 0, b2_T(i)));
                                break
                            end
                        end
                    end
                end
                
            elseif j == 2 || j == 4
                if i == 1
                    if j == 4
                        e2_T_final(i) = double(solve(A2_eq{4}(1) == 0));
                    end
                elseif i > 1
                    if b2_T_final(i-1) > 1e-7
                        for k = 1:2:A2_eq_column(j)
                            if double(solve(A2_eq{j}(k) == 0, e2_T(i))) > 1e-7
                                e2_T_final(i) = double(solve(A2_eq{j}(k) == 0, e2_T(i)));
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    A2_matrix{i} = [1 b2_T_final(i)*s; e2_T_final(i)*s 1];
    T2_temp = [1 -b2_T_final(i)*s; -e2_T_final(i)*s 1]*T2_temp;
    if i < 2*N
        A2_inv_T = [1 -b2_T(i+1)*s; -e2_T(i+1)*s 1]*T2_temp;
    end
end

b2_T_final
e2_T_final

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Denormalizing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_denorm = RL_denorm;

L1_vals = b1_T_final*K_denorm
C1_vals = e1_T_final/K_denorm

L2_vals = b2_T_final*K_denorm
C2_vals = e2_T_final/K_denorm

if Rg/RL == 1
    L1_vals_dual = e1_T_final*K_denorm
    C1_vals_dual = b1_T_final/K_denorm

    L2_vals_dual = e2_T_final*K_denorm
    C2_vals_dual = b2_T_final/K_denorm
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% second method to calculate p for F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RL2 = RL_denorm/Rg_denorm;
Rg2 = 1;

[DEN1_columns, DEN1_rows] = size(DEN1);
K1 = DEN1(DEN1_rows)*RL2/(RL2+Rg2);
DEN_F1_s = poly2sym(DEN1, s);
p1_squared = 1 - 4/RL2*K1/DEN_F1_s*K1/subs(DEN_F1_s, s, -s);
p1_squared_NUM = sym2poly(simplify(p1_squared*DEN_F1_s*subs(DEN_F1_s,s,-s)));
p1_squared_zeros = roots(p1_squared_NUM);

sum = 0;
for i = 1:size(p1_squared_zeros,1)
    sum = sum + p1_squared_zeros(i);
end

if sum > 1e-7
    fprintf("F1 has no more LC circuits \r");
else
    p1_zeros2 = p1_squared_zeros(abs(real(p1_squared_zeros))<1e-10);
    p1_zeros3 = p1_zeros2(abs(p1_zeros2)<1e-10);
    if mod(size(p1_zeros2, 1),4) ~= 0 && mod(size(p1_zeros3, 1),2) ~= 0
        fprintf("F1 has no more LC circuits \r");
    else
        if size(p1_zeros2, 1) ~= 0
            p1_zeros5 = zeros(1, 2*N);
            for j = 1:size(p1_zeros2, 1)/4
                for i = 1:size(p1_zeros2, 1)
                    if abs(real(p1_zeros2(i)))<1e-10
                        p1_zeros5(i+2*(j-1)) = p1_zeros2(1+4*(j-1));
                        p1_zeros5(i+1+2*(j-1)) = p1_zeros2(2+4*(j-1));
                    end
                end
            end
        end

        p1_zeros4 = p1_squared_zeros(real(p1_squared_zeros)<0);
        j = 1;
        for i = size(p1_zeros2, 1)/4*2+1:2*N
            p1_zeros5(i) = p1_zeros4(j);
            j = j + 1;
        end

        [p1_zeros5_rows, p1_zeros5_columns] = size(p1_zeros5);
        for i = 1:p1_zeros5_columns
            if abs(real(p1_zeros5(i)))<1e-7
                p1_zeros5(i) = imag(p1_zeros5(i))*1i;
            elseif abs(imag(p1_zeros5(i)))<1e-7
                p1_zeros5(i) = real(p1_zeros5(i));
            elseif abs(real(p1_zeros5(i)))<1e-7 && abs(imag(p1_zeros5(i)))<1e-7
                p1_zeros5(i) = 0;
            end
        end

        p1_NUM = [1 -p1_zeros5(1)];
        for i = 2:2*N
            p1_NUM_temp = [1 -p1_zeros5(i)];
            p1_NUM = conv(p1_NUM, p1_NUM_temp);
        end

        [p1_NUM_rows, p1_NUM_columns] = size(p1_NUM);
        for i = 1:p1_NUM_columns
            if abs(imag(p1_NUM(i)))<1e-7
                p1_NUM(i) = real(p1_NUM(i));
            elseif abs(real(p1_NUM(i)))<1e-7 && abs(imag(p1_NUM(i)))<1e-7
                p1_NUM(i) = 0;
            end
        end

        p1 = vpa(simplify(poly2sym(p1_NUM, s)/poly2sym(DEN_F1_s, s)));
    end
end

if RL2==1
    fprintf("F1 has dual \r");
    Zin1 = vpa(simplify((1+p1)/(1-p1)));
else
    if eval(subs((1+p1)/(1-p1), s, 0)) == RL2
        Zin1 = vpa(simplify((1+p1)/(1-p1)));
    elseif eval(subs((1+p1)/(1-p1), s, 0)) == 1/RL2
        Zin1 = vpa(simplify((1-p1)/(1+p1)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cauer LC Zin1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_Zin1,den_Zin1] = numden(Zin1);
num_Zin1 = sym2poly(num_Zin1);
den_Zin1 = sym2poly(den_Zin1);
tf(num_Zin1,den_Zin1)
Zin1_func_is_suitable_for_cauer1 = true;
for i=1:length(num_Zin1)
    if num_Zin1(i) ~= real(num_Zin1(i))
        Zin1_func_is_suitable_for_cauer1 = false;
        disp("cauer1 can not be performed with Zin1")
    end
end
if Zin1_func_is_suitable_for_cauer1
    for i=1:length(den_Zin1)
        if den_Zin1(i) ~= real(den_Zin1(i))
            Zin1_func_is_suitable_for_cauer1 = false;
            disp("cauer1 can not be performed with Zin1")
            break
        end
    end
end
if Zin1_func_is_suitable_for_cauer1
    %perform cauer 1 on Zin1
    kn_Zin1 = zeros(1,20);
    new_num_FLCs_Zin1 = num_Zin1;
    new_den_FLCs_Zin1 = den_Zin1;
    new_den_sFLCs_Zin1 = conv([1 0],new_den_FLCs_Zin1); %sFLCs stand for (1/s)*FLCs
    sys_syms_Zin1 = poly2sym(new_num_FLCs_Zin1,s)/poly2sym(new_den_sFLCs_Zin1,s); % limit sFLCs in s->inf
    kn_Zin1(1) = limit(sys_syms_Zin1, s, inf);
    if kn_Zin1(1)<0
        disp("cauer1 can not be performed with Zin1")
        Zin1_func_is_suitable_for_cauer1 = false;
    end
    num_FLCs = new_den_FLCs_Zin1;
    den_FLCs = new_num_FLCs_Zin1 - (kn_Zin1(1)*new_den_sFLCs_Zin1);
    den_sFLCs = conv([1 0], den_FLCs);
    den_FLCs(den_FLCs<1e-2) = 0;
    den_FLCs
    %delete zeros at the begining of den_FLCs and den_sFLCs
    counter1 = 0;
    counter2 = 0;
    for j=1:length(den_FLCs)
        if den_FLCs(j)==0
            counter1 = counter1 + 1;
        else
            break
        end
    end

    for jj=1:length(den_sFLCs)
        if den_sFLCs(jj)==0
            counter2 = counter2 + 1;
        else
            break
        end
    end           
    den_FLCs = den_FLCs(counter1+1:end);
    den_sFLCs = den_sFLCs(counter2+1:end);
end

%while
i = 2;
   while (length(den_FLCs)~=1) && Zin1_func_is_suitable_for_cauer1
       if (length(num_FLCs) == length(den_FLCs)) && (rem(i,2)==0) %deconv process
           %disp('deconv started...')
           [q,r] = deconv(num_FLCs,den_FLCs);
           kn_Zin1(i) = q;
           if kn_Zin1(i)<0
               disp("cauer1 can not be performed with Zin1")
               Zin1_func_is_suitable_for_cauer1 = false;
           end
           num_FLCs = den_FLCs;
           den_FLCs = r;
           den_sFLCs = conv([1 0],den_FLCs);
           % delete zeros in the begining of den_FLCs and den_sFLCs
           counter1 = 0;
           counter2 = 0;
           for j=1:length(den_FLCs)
                if den_FLCs(j)==0
                    counter1 = counter1 + 1;
                else
                    break
                end
            end

            for jj=1:length(den_sFLCs)
                if den_sFLCs(jj)==0
                    counter2 = counter2 + 1;
                else
                    break
                end
            end
            den_FLCs = den_FLCs(counter1+1:end);
            den_sFLCs = den_sFLCs(counter2+1:end);
            
          
       elseif (length(num_FLCs) == length(den_FLCs)+1) && (rem(i,2)==1)   %cauer1 process
            %disp("cauer1 preocess started...")
                new_num_FLCs_Zin1 = num_FLCs;
                new_den_FLCs_Zin1 = den_FLCs;
                new_den_sFLCs_Zin1 = den_sFLCs;
                
                sys_syms_Zin1 = poly2sym(new_num_FLCs_Zin1,s)/poly2sym(new_den_sFLCs_Zin1,s);
                kn_Zin1(i) = limit(sys_syms_Zin1, s, inf);
                if kn_Zin1(i)<0
                    disp("cauer1 can not be performed with Zin1")
                    Zin1_func_is_suitable_for_cauer1 = false;
                end
                
                num_FLCs = new_den_FLCs_Zin1;
                den_FLCs = new_num_FLCs_Zin1 - (kn_Zin1(i).*new_den_sFLCs_Zin1);
                den_sFLCs = conv([1 0],den_FLCs);
                
                num_FLCs(abs(num_FLCs)<0.000001) = 0; %khataye matlab
                den_FLCs(abs(den_FLCs)<0.000001) = 0;
                den_sFLCs(abs(den_sFLCs)<0.000001) = 0;
                
                %delete zeros at the begining of den_FLCs and den_sFLCs
                counter1 = 0;
                counter2 = 0;
                for j=1:length(den_FLCs)
                    if den_FLCs(j)==0
                        counter1 = counter1 + 1;
                    else
                        break
                    end
                end

                for jj=1:length(den_sFLCs)
                    if den_sFLCs(jj)==0
                        counter2 = counter2 + 1;
                    else
                        break
                    end
                end
           
                den_FLCs = den_FLCs(counter1+1:end);
                den_sFLCs = den_sFLCs(counter2+1:end);

       end
       i = i+1;
   end
   if Zin1_func_is_suitable_for_cauer1
        kn_Zin1(i) = num_FLCs(1)/den_FLCs;
   end
    if kn_Zin1(i)<0
        disp("cauer1 can not be performed with Zin1")
        Zin1_func_is_suitable_for_cauer1 = false;
    end
   if Zin1_func_is_suitable_for_cauer1
        i = i+1;
        kn_Zin1(i) = num_FLCs(2)/den_FLCs;
   end
       if kn_Zin1(i)<0
            disp("cauer1 can not be performed with Zin1")
            Zin1_func_is_suitable_for_cauer1 = false;
       end
   %hazf sefrhaye baad tahe kn_F1
   for iii=1:length(kn_Zin1)
       if kn_Zin1(iii) == 0
           break
       end
   end
   kn_Zin1 = kn_Zin1(1:iii-1);
 if Zin1_func_is_suitable_for_cauer1
        kn_Zin1
 end

   
   % seperation of kn_Zin1 to L and C
   C_y22_Zin1_kuh = zeros(1,floor(length(kn_Zin1)/2));
   L_y22_Zin1_kuh = zeros(1,floor(length(kn_Zin1)/2));
   index_L = 1;
   index_C =1;
   for iter=1:length(kn_Zin1)
       if rem(iter,2) == 1
           L_y22_Zin1_kuh(index_L) = kn_Zin1(iter);
           index_L = index_L + 1;
       elseif rem(iter,2) == 0
           C_y22_Zin1_kuh(index_C) = kn_Zin1(iter);
           index_C = index_C + 1;
       end
   end
 if Zin1_func_is_suitable_for_cauer1
        kn_Zin1
        C_y22_Zin1_kuh
        L_y22_Zin1_kuh
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% second method to calculate p for F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DEN2_columns, DEN2_rows] = size(DEN2);
K2 = DEN2(DEN2_rows)*RL2/(RL2+Rg2);
DEN_F2_s = poly2sym(DEN2, s);
p2_squared = 1 - 4/RL2*K2/DEN_F2_s*K2/subs(DEN_F2_s, s, -s);
p2_squared_NUM = sym2poly(simplify(p2_squared*DEN_F2_s*subs(DEN_F2_s,s,-s)));
p2_squared_zeros = roots(p2_squared_NUM);

sum = 0;
for i = 1:size(p2_squared_zeros,1)
    sum = sum + p2_squared_zeros(i);
end

if sum > 1e-7
    fprintf("F1 has no more LC circuits \r");
else
    p2_zeros2 = p2_squared_zeros(abs(real(p2_squared_zeros))<1e-7);
    p2_zeros3 = p2_zeros2(abs(p2_zeros2)<1e-10);

    p2_zeros5 = zeros(1, 2*N);
    p2_zeros6 = sort(p2_zeros2);
    j = 1;
    for i = 0:4:size(p2_zeros6,1)-1
        p2_zeros5(j) = p2_zeros6(i+1);
        p2_zeros5(j+1) = p2_zeros6(i+2);
        j = j+2;
    end

    [p2_zeros5_rows, p2_zeros5_columns] = size(p2_zeros5);
    for i = 1:p2_zeros5_columns
        if abs(real(p2_zeros5(i)))<1e-7
            p2_zeros5(i) = imag(p2_zeros5(i))*1i;
        elseif abs(imag(p2_zeros5(i)))<1e-7
            p2_zeros5(i) = real(p2_zeros5(i));
        elseif abs(real(p2_zeros5(i)))<1e-7 && abs(imag(p2_zeros5(i)))<1e-7
            p2_zeros5(i) = 0;
        end
    end

    p2_NUM = [1 -p2_zeros5(1)];
    for i = 2:2*N
        p2_NUM_temp = [1 -p2_zeros5(i)];
        p2_NUM = conv(p2_NUM, p2_NUM_temp);
    end

    [p2_NUM_rows, p2_NUM_columns] = size(p2_NUM);
    for i = 1:p2_NUM_columns
        if abs(imag(p2_NUM(i)))<1e-7
            p2_NUM(i) = real(p2_NUM(i));
        elseif abs(real(p2_NUM(i)))<1e-7 && abs(imag(p2_NUM(i)))<1e-7
            p2_NUM(i) = 0;
        end
    end

    p2 = vpa(simplify(poly2sym(p2_NUM, s)/poly2sym(DEN_F2_s, s)));

    if RL2==1
        fprintf("F2 has dual \r");
        Zin2 = vpa(simplify((1+p2)/(1-p2)));
    else
        if eval(subs((1+p2)/(1-p2), s, 0)) == RL2
            Zin2 = vpa(simplify((1+p2)/(1-p2)));
        elseif eval(subs((1+p2)/(1-p2), s, 0)) == 1/RL2
            Zin2 = vpa(simplify((1-p2)/(1+p2)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cauer LC Zin2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_Zin2,den_Zin2] = numden(Zin2);
num_Zin2 = sym2poly(num_Zin2);
den_Zin2 = sym2poly(den_Zin2);

tf(num_Zin2,den_Zin2)
Zin2_func_is_suitable_for_cauer1 = true;
for i=1:length(num_Zin2)
    if num_Zin2(i) ~= real(num_Zin2(i))
        Zin2_func_is_suitable_for_cauer1 = false;
        disp("cauer1 can not be performed with Zin2")
        break
    end
end
if Zin2_func_is_suitable_for_cauer1
    for i=1:length(den_Zin2)
        if den_Zin2(i) ~= real(den_Zin2(i))
            Zin2_func_is_suitable_for_cauer1 = false;
            disp("cauer1 can not be performed with Zin2")
            break
        end
    end
end
if Zin2_func_is_suitable_for_cauer1
    %perform cauer 1 on Zin2
    kn_Zin2 = zeros(1,20);
    new_num_FLCs_Zin2 = num_Zin2;
    new_den_FLCs_Zin2 = den_Zin2;
    new_den_sFLCs_Zin2 = conv([1 0],new_den_FLCs_Zin2); %sFLCs stand for (1/s)*FLCs
    sys_syms_Zin2 = poly2sym(new_num_FLCs_Zin2,s)/poly2sym(new_den_sFLCs_Zin2,s); % limit sFLCs in s->inf
    kn_Zin2(1) = limit(sys_syms_Zin2, s, inf);
    if kn_Zin2(1)<0
        disp("cauer1 can not be performed with Zin2")
        Zin2_func_is_suitable_for_cauer1 = false;
    end
    num_FLCs = new_den_FLCs_Zin2;
    den_FLCs = new_num_FLCs_Zin2 - (kn_Zin2(1)*new_den_sFLCs_Zin2);
    den_sFLCs = conv([1 0], den_FLCs);
    den_FLCs(abs(den_FLCs)<1e-2) = 0;

    %delete zeros at the begining of den_FLCs and den_sFLCs
    counter1 = 0;
    counter2 = 0;
    for j=1:length(den_FLCs)
        if den_FLCs(j)==0
            counter1 = counter1 + 1;
        else
            break
        end
    end

    for jj=1:length(den_sFLCs)
        if den_sFLCs(jj)==0
            counter2 = counter2 + 1;
        else
            break
        end
    end           
    den_FLCs = den_FLCs(counter1+1:end);
    den_sFLCs = den_sFLCs(counter2+1:end);
end

%while
i = 2;
   while (length(den_FLCs)~=1) && Zin2_func_is_suitable_for_cauer1
       if (length(num_FLCs) == length(den_FLCs)) && (rem(i,2)==0) %deconv process
           %disp('deconv started...')
           [q,r] = deconv(num_FLCs,den_FLCs);
           kn_Zin2(i) = q;
           if kn_Zin2(i)<0
               disp("cauer1 can not be performed with Zin2")
               Zin2_func_is_suitable_for_cauer1 = false;
           end
           num_FLCs = den_FLCs;
           den_FLCs = r;
           den_sFLCs = conv([1 0],den_FLCs);
           % delete zeros in the begining of den_FLCs and den_sFLCs
           counter1 = 0;
           counter2 = 0;
           for j=1:length(den_FLCs)
                if den_FLCs(j)==0
                    counter1 = counter1 + 1;
                else
                    break
                end
            end

            for jj=1:length(den_sFLCs)
                if den_sFLCs(jj)==0
                    counter2 = counter2 + 1;
                else
                    break
                end
            end
            den_FLCs = den_FLCs(counter1+1:end);
            den_sFLCs = den_sFLCs(counter2+1:end);
            
          
       elseif (length(num_FLCs) == length(den_FLCs)+1) && (rem(i,2)==1)   %cauer1 process
            %disp("cauer1 preocess started...")
                new_num_FLCs_Zin2 = num_FLCs;
                new_den_FLCs_Zin2 = den_FLCs;
                new_den_sFLCs_Zin2 = den_sFLCs;
                
                sys_syms_Zin2 = poly2sym(new_num_FLCs_Zin2,s)/poly2sym(new_den_sFLCs_Zin2,s);
                kn_Zin2(i) = limit(sys_syms_Zin2, s, inf);
                if kn_Zin2(i)<0
                    disp("cauer1 can not be performed with Zin2")
                    Zin2_func_is_suitable_for_cauer1 = false;
                end
                
                num_FLCs = new_den_FLCs_Zin2;
                den_FLCs = new_num_FLCs_Zin2 - (kn_Zin2(i).*new_den_sFLCs_Zin2);
                den_sFLCs = conv([1 0],den_FLCs);
                
                num_FLCs(abs(num_FLCs)<0.000001) = 0; %khataye matlab
                den_FLCs(abs(den_FLCs)<0.000001) = 0;
                den_sFLCs(abs(den_sFLCs)<0.000001) = 0;
                
                %delete zeros at the begining of den_FLCs and den_sFLCs
                counter1 = 0;
                counter2 = 0;
                for j=1:length(den_FLCs)
                    if den_FLCs(j)==0
                        counter1 = counter1 + 1;
                    else
                        break
                    end
                end

                for jj=1:length(den_sFLCs)
                    if den_sFLCs(jj)==0
                        counter2 = counter2 + 1;
                    else
                        break
                    end
                end
           
                den_FLCs = den_FLCs(counter1+1:end);
                den_sFLCs = den_sFLCs(counter2+1:end);

       end
       i = i+1;
   end
   if Zin2_func_is_suitable_for_cauer1
        kn_Zin2(i) = num_FLCs(1)/den_FLCs;
        if kn_Zin2(i)<0
            disp("cauer1 can not be performed with Zin2")
            Zin2_func_is_suitable_for_cauer1 = false;
        end
   end

   if Zin2_func_is_suitable_for_cauer1
        i = i+1;
        kn_Zin2(i) = num_FLCs(2)/den_FLCs;
        if kn_Zin2(i)<0
            disp("cauer1 can not be performed with Zin2")
            Zin2_func_is_suitable_for_cauer1 = false;
        end
   end

   %hazf sefrhaye baad tahe kn_Zin2
   if Zin2_func_is_suitable_for_cauer1
       for iii=1:length(kn_Zin2)
           if kn_Zin2(iii) == 0
               break
           end
       end
       kn_Zin2 = kn_Zin2(1:iii-1);
     if Zin2_func_is_suitable_for_cauer1
            kn_Zin2
     end


       % seperation of kn_Zin2 to L and C
       C_y22_Zin2_kuh = zeros(1,floor(length(kn_Zin2)/2));
       L_y22_Zin2_kuh = zeros(1,floor(length(kn_Zin2)/2));
       index_L = 1;
       index_C =1;
       for iter=1:length(kn_Zin2)
           if rem(iter,2) == 1
               L_y22_Zin2_kuh(index_L) = kn_Zin2(iter);
               index_L = index_L + 1;
           elseif rem(iter,2) == 0
               C_y22_Zin2_kuh(index_C) = kn_Zin2(iter);
               index_C = index_C + 1;
           end
       end
     if Zin2_func_is_suitable_for_cauer1
            kn_Zin2
            C_y22_Zin2_kuh
            L_y22_Zin2_kuh
     end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yanagisawa F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nD1_y = 2*N-1;

D1_y = [1 0];
for i = 1:nD1_y
    D1_y_temp = [1 i];
    D1_y = conv(D1_y, D1_y_temp);
end

[r_A1_sD1, p_A1_sD1, k_A1_sD1] = residue(NUM1, D1_y);
[r_B1_sD1, p_B1_sD1, k_B1_sD1] = residue(DEN1, D1_y);

ya1_y = 0;
yb1_y = 0;
for i = 1:size(r_A1_sD1, 1)
    if r_A1_sD1(i)>0
        ya1_y = ya1_y + r_A1_sD1(i)/(s-p_A1_sD1(i));
    elseif r_A1_sD1(i)<0
        yb1_y = yb1_y - r_A1_sD1(i)/(s-p_A1_sD1(i));
    end
end

if size(k_A1_sD1, 1) ~= 0
    for i = 1:size(k_A1_sD1, 1)
        if k_A1_sD1(i)>0
            ya1_y = ya1_y + k_A1_sD1(i)*s*(i-1);
        elseif k_A1_sD1(i)<0
            yb1_y = yb1_y - k_A1_sD1(i)*s*(i-1);
        end
    end
end

yc1_y = 0;
ya1_yb1_yd1 = 0;
for i = 1:size(r_B1_sD1, 1)
    if r_B1_sD1(i)>0
        yc1_y = yc1_y + r_B1_sD1(i)/(s-p_B1_sD1(i));
    elseif r_B1_sD1(i)<0
        ya1_yb1_yd1 = ya1_yb1_yd1 - r_B1_sD1(i)/(s-p_B1_sD1(i));
    end
end

if size(k_B1_sD1, 1) ~= 0
    for i = 1:size(k_B1_sD1, 1)
        if k_B1_sD1(i)>0
            yc1_y = yc1_y + k_B1_sD1(i)*s*(i-1);
        elseif k_B1_sD1(i)<0
            ya1_yb1_yd1 = ya1_yb1_yd1 - k_B1_sD1(i)*s*(i-1);
        end
    end
end
yd1_y = ya1_y + yb1_y - ya1_yb1_yd1;

ya1_y = vpa(simplify(s*ya1_y));
K11_y = 1;
yb1_y = vpa(simplify(s*yb1_y));
yc1_y = vpa(simplify(s*yc1_y));
K21_y = 1;
yd1_y = vpa(simplify(s*yd1_y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yanagisawa F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nD2_y = 2*N-1;

D2_y = [1 0];
for i = 1:nD2_y
    D2_y_temp = [1 i];
    D2_y = conv(D2_y, D2_y_temp);
end

[r_A2_sD2, p_A2_sD2, k_A2_sD2] = residue(NUM2, D2_y);
[r_B2_sD2, p_B2_sD2, k_B2_sD2] = residue(DEN2, D2_y);

ya2_y = 0;
yb2_y = 0;
for i = 1:size(r_A2_sD2, 1)
    if r_A2_sD2(i)>0
        ya2_y = ya2_y + r_A2_sD2(i)/(s-p_A2_sD2(i));
    elseif r_A2_sD2(i)<0
        yb2_y = yb2_y - r_A2_sD2(i)/(s-p_A2_sD2(i));
    end
end

if size(k_A2_sD2, 1) ~= 0
    for i = 1:size(k_A2_sD2, 1)
        if k_A2_sD2(i)>0
            ya2_y = ya2_y + k_A2_sD2(i)*s*(i-1);
        elseif k_A2_sD2(i)<0
            yb2_y = yb2_y - k_A2_sD2(i)*s*(i-1);
        end
    end
end

yc2_y = 0;
ya1_yb1_yd1 = 0;
for i = 1:size(r_B2_sD2, 1)
    if r_B2_sD2(i)>0
        yc2_y = yc2_y + r_B2_sD2(i)/(s-p_B2_sD2(i));
    elseif r_B2_sD2(i)<0
        ya1_yb1_yd1 = ya1_yb1_yd1 + r_B2_sD2(i)/(s-p_B2_sD2(i));
    end
end

if size(k_B2_sD2, 1) ~= 0
    for i = 1:size(k_B2_sD2, 1)
        if k_B2_sD2(i)>0
            yc2_y = yc2_y + k_B2_sD2(i)*s*(i-1);
        elseif k_B2_sD2(i)<0
            ya1_yb1_yd1 = ya1_yb1_yd1 + k_B2_sD2(i)*s*(i-1);
        end
    end
end
yd2_y = ya2_y + yb2_y - ya1_yb1_yd1;

ya2_y = vpa(simplify(s*ya2_y));
K12_y = 1;
yb2_y = vpa(simplify(s*yb2_y));
yc2_y = vpa(simplify(s*yc2_y));
K22_y = 1;
yd2_y = vpa(simplify(s*yd2_y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mathews Seifret F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ya1_ms_1 = 0;
yb1_ms_1 = 0;
ya1_ms_2 = 0;
yb1_ms_2 = 0;
for i = 1:size(r_A1_sD1, 1)
    if r_A1_sD1(i)>0
        ya1_ms_1 = ya1_ms_1 + r_A1_sD1(i)/(s-p_A1_sD1(i));
        yb1_ms_2 = yb1_ms_2 + r_A1_sD1(i)/(s-p_A1_sD1(i));
    elseif r_A1_sD1(i)<0
        ya1_ms_2 = ya1_ms_2 - r_A1_sD1(i)/(s-p_A1_sD1(i));
        yb1_ms_1 = yb1_ms_1 - r_A1_sD1(i)/(s-p_A1_sD1(i));
    end
end

if size(k_A1_sD1, 1) ~= 0
    for i = 1:size(k_A1_sD1, 1)
        if k_A1_sD1(i)>0
            ya1_ms_1 = ya1_ms_1 + k_A1_sD1(i)*s*(i-1);
            yb1_ms_2 = yb1_ms_2 + k_A1_sD1(i)*s*(i-1);
        elseif k_A1_sD1(i)<0
            ya1_ms_2 = ya1_ms_2 - k_A1_sD1(i)*s*(i-1);
            yb1_ms_1 = yb1_ms_1 - k_A1_sD1(i)*s*(i-1);
        end
    end
end

yc1_ms_1 = 0;
yd1_ms_1 = 0;
yc1_ms_2 = 0;
yd1_ms_2 = 0;
for i = 1:size(r_B1_sD1, 1)
    if r_B1_sD1(i)>0
        yc1_ms_2 = yc1_ms_2 + r_B1_sD1(i)/(s-p_B1_sD1(i));
        yd1_ms_1 = yd1_ms_1 + r_B1_sD1(i)/(s-p_B1_sD1(i));
    elseif r_B1_sD1(i)<0
        yc1_ms_1 = yc1_ms_1 - r_B1_sD1(i)/(s-p_B1_sD1(i));
        yd1_ms_2 = yd1_ms_2 - r_B1_sD1(i)/(s-p_B1_sD1(i));
    end
end

if size(k_B1_sD1, 1) ~= 0
    for i = 1:size(k_B1_sD1, 1)
        if k_B1_sD1(i)>0
            yc1_ms_2 = yc1_ms_2 + k_B1_sD1*s*(i-1);
            yd1_ms_1 = yd1_ms_1 + k_B1_sD1*s*(i-1);
        elseif k_B1_sD1(i)<0
            yc1_ms_1 = yc1_ms_1 - k_B1_sD1*s*(i-1);
            yd1_ms_2 = yd1_ms_2 - k_B1_sD1*s*(i-1);
        end
    end
end

if yc1_ms_1 == 0
    yc1_ms = vpa(simplify(s*yc1_ms_2));
    yd1_ms = vpa(simplify(s*yd1_ms_2));
    ya1_ms = vpa(simplify(s*ya1_ms_2));
    yb1_ms = vpa(simplify(s*yb1_ms_2));
elseif yc1_ms_2 == 0
    yc1_ms = vpa(simplify(s*yc1_ms_1));
    yd1_ms = vpa(simplify(s*yd1_ms_1));
    ya1_ms = vpa(simplify(s*ya1_ms_1));
    yb1_ms = vpa(simplify(s*yb1_ms_1));
else
    yc1_ms = vpa(simplify(s*yc1_ms_1));
    yd1_ms = vpa(simplify(s*yd1_ms_1));
    ya1_ms = vpa(simplify(s*ya1_ms_1));
    yb1_ms = vpa(simplify(s*yb1_ms_1));
end

K11_ms = 1;
K21_ms = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mathews Seifret F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ya2_ms_1 = 0;
yb2_ms_1 = 0;
ya2_ms_2 = 0;
yb2_ms_2 = 0;
for i = 1:size(r_A2_sD2, 1)
    if r_A2_sD2(i)>0
        ya2_ms_1 = ya2_ms_1 + r_A2_sD2(i)/(s-p_A2_sD2(i));
        yb2_ms_2 = yb2_ms_2 + r_A2_sD2(i)/(s-p_A2_sD2(i));
    elseif r_A2_sD2(i)<0
        ya2_ms_2 = ya2_ms_2 - r_A2_sD2(i)/(s-p_A2_sD2(i));
        yb2_ms_1 = yb2_ms_1 - r_A2_sD2(i)/(s-p_A2_sD2(i));
    end
end

if size(k_A2_sD2, 1) ~= 0
    for i = 1:size(k_A2_sD2, 1)
        if k_A2_sD2(i)>0
            ya2_ms_1 = ya2_ms_1 + k_A2_sD2(i)*s*(i-1);
            yb2_ms_2 = yb2_ms_2 + k_A2_sD2(i)*s*(i-1);
        elseif k_A2_sD2(i)<0
            ya2_ms_2 = ya2_ms_2 - k_A2_sD2(i)*s*(i-1);
            yb2_ms_1 = yb2_ms_1 - k_A2_sD2(i)*s*(i-1);
        end
    end
end

yc2_ms_1 = 0;
yd2_ms_1 = 0;
yc2_ms_2 = 0;
yd2_ms_2 = 0;
for i = 1:size(r_B2_sD2, 1)
    if r_B2_sD2(i)>0
        yc2_ms_2 = yc2_ms_2 + r_B2_sD2(i)/(s-p_B2_sD2(i));
        yd2_ms_1 = yd2_ms_1 + r_B2_sD2(i)/(s-p_B2_sD2(i));
    elseif r_B2_sD2(i)<0
        yc2_ms_1 = yc2_ms_1 - r_B2_sD2(i)/(s-p_B2_sD2(i));
        yd2_ms_2 = yd2_ms_2 - r_B2_sD2(i)/(s-p_B2_sD2(i));
    end
end

if size(k_B2_sD2, 1) ~= 0
    for i = 1:size(k_B2_sD2, 1)
        if k_B2_sD2(i)>0
            yc2_ms_2 = yc2_ms_2 + k_B2_sD2*s*(i-1);
            yd2_ms_1 = yd2_ms_1 + k_B2_sD2*s*(i-1);
        elseif k_B2_sD2(i)<0
            yc2_ms_1 = yc2_ms_1 - k_B2_sD2*s*(i-1);
            yd2_ms_2 = yd2_ms_2 - k_B2_sD2*s*(i-1);
        end
    end
end

if yc2_ms_1 == 0
    yc2_ms = vpa(simplify(s*yc2_ms_2));
    yd2_ms = vpa(simplify(s*yd2_ms_2));
    ya2_ms = vpa(simplify(s*ya2_ms_2));
    yb2_ms = vpa(simplify(s*yb2_ms_2));
elseif yc2_ms_2 == 0
    yc2_ms = vpa(simplify(s*yc2_ms_1));
    yd2_ms = vpa(simplify(s*yd2_ms_1));
    ya2_ms = vpa(simplify(s*ya2_ms_1));
    yb2_ms = vpa(simplify(s*yb2_ms_1));
else
    yc2_ms = vpa(simplify(s*yc2_ms_1));
    yd2_ms = vpa(simplify(s*yd2_ms_1));
    ya2_ms = vpa(simplify(s*ya2_ms_1));
    yb2_ms = vpa(simplify(s*yb2_ms_1));
end

K12_ms = 1;
K22_ms = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lovering F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yc1_l = yc1_ms;
yd1_l = yd1_ms;
ya1_1 = ya1_ms;
yb1_l = yb1_ms;

yo1_l = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lovering F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yc2_l = yc2_ms;
yd2_l = yd2_ms;
ya2_1 = ya2_ms;
yb2_l = yb2_ms;

yo2_l = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mitra F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ya1_m_1 = 0;
yb1_m_1 = 0;
ya1_m_2 = 0;
yb1_m_2 = 0;
for i = 1:size(r_A1_sD1, 1)
    if r_A1_sD1(i)>0
        ya1_m_1 = ya1_m_1 + r_A1_sD1(i)/(s-p_A1_sD1(i));
        yb1_m_2 = yb1_m_2 + r_A1_sD1(i)/(s-p_A1_sD1(i));
    elseif r_A1_sD1(i)<0
        ya1_m_2 = ya1_m_2 - r_A1_sD1(i)/(s-p_A1_sD1(i));
        yb1_m_1 = yb1_m_1 - r_A1_sD1(i)/(s-p_A1_sD1(i));
    end
end

if size(k_A1_sD1, 1) ~= 0
    for i = 1:size(k_A1_sD1, 1)
        if k_A1_sD1(i)>0
            ya1_m_1 = ya1_m_1 + k_A1_sD1(i)*s*(i-1);
            yb1_m_2 = yb1_m_2 + k_A1_sD1(i)*s*(i-1);
        elseif k_A1_sD1(i)<0
            ya1_m_2 = ya1_m_2 - k_A1_sD1(i)*s*(i-1);
            yb1_m_1 = yb1_m_1 - k_A1_sD1(i)*s*(i-1);
        end
    end
end

ye1_m_1 = 0;
yf1_m_1 = 0;
ye1_m_2 = 0;
yf1_m_2 = 0;
for i = 1:size(r_B1_sD1, 1)
    if r_B1_sD1(i)>0
        ye1_m_2 = ye1_m_2 + r_B1_sD1(i)/(s-p_B1_sD1(i));
        yf1_m_1 = yf1_m_1 + r_B1_sD1(i)/(s-p_B1_sD1(i));
    elseif r_B1_sD1(i)<0
        ye1_m_1 = ye1_m_1 - r_B1_sD1(i)/(s-p_B1_sD1(i));
        yf1_m_2 = yf1_m_2 - r_B1_sD1(i)/(s-p_B1_sD1(i));
    end
end

if size(k_B1_sD1, 1) ~= 0
    for i = 1:size(k_B1_sD1, 1)
        if k_B1_sD1(i)>0
            ye1_m_2 = ye1_m_2 + k_B1_sD1*s*(i-1);
            yf1_m_1 = yf1_m_1 + k_B1_sD1*s*(i-1);
        elseif k_B1_sD1(i)<0
            ye1_m_1 = ye1_m_1 - k_B1_sD1*s*(i-1);
            yf1_m_2 = yf1_m_2 - k_B1_sD1*s*(i-1);
        end
    end
end

if yf1_m_1 == 0
    yf1_m = vpa(simplify(s*ye1_m_2));
    ye1_m = vpa(simplify(s*yf1_m_2));
    ya1_m = vpa(simplify(s*ya1_m_2));
    yb1_m = vpa(simplify(s*yb1_m_2));
elseif yf1_m_2 == 0
    yf1_m = vpa(simplify(s*ye1_m_1));
    ye1_m = vpa(simplify(s*yf1_m_1));
    ya1_m = vpa(simplify(s*ya1_m_1));
    yb1_m = vpa(simplify(s*yb1_m_1));
else
    yf1_m = vpa(simplify(s*ye1_m_1));
    ye1_m = vpa(simplify(s*yf1_m_1));
    ya1_m = vpa(simplify(s*ya1_m_1));
    yb1_m = vpa(simplify(s*yb1_m_1));
end

yc1_m = yb1_m + yf1_m;
yd1_m = ya1_m + ye1_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mitra F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ya2_m_1 = 0;
yb2_m_1 = 0;
ya2_m_2 = 0;
yb2_m_2 = 0;
for i = 1:size(r_A2_sD2, 1)
    if r_A2_sD2(i)>0
        ya2_m_1 = ya2_m_1 + r_A2_sD2(i)/(s-p_A2_sD2(i));
        yb2_m_2 = yb2_m_2 + r_A2_sD2(i)/(s-p_A2_sD2(i));
    elseif r_A2_sD2(i)<0
        ya2_m_2 = ya2_m_2 - r_A2_sD2(i)/(s-p_A2_sD2(i));
        yb2_m_1 = yb2_m_1 - r_A2_sD2(i)/(s-p_A2_sD2(i));
    end
end

if size(k_A2_sD2, 1) ~= 0
    for i = 1:size(k_A2_sD2, 1)
        if k_A2_sD2(i)>0
            ya2_m_1 = ya2_m_1 + k_A2_sD2(i)*s*(i-1);
            yb2_m_2 = yb2_m_2 + k_A2_sD2(i)*s*(i-1);
        elseif k_A2_sD2(i)<0
            ya2_m_2 = ya2_m_2 - k_A2_sD2(i)*s*(i-1);
            yb2_m_1 = yb2_m_1 - k_A2_sD2(i)*s*(i-1);
        end
    end
end

ye2_m_1 = 0;
yf2_m_1 = 0;
ye2_m_2 = 0;
yf2_m_2 = 0;
for i = 1:size(r_B2_sD2, 1)
    if r_B2_sD2(i)>0
        ye2_m_2 = ye2_m_2 + r_B2_sD2(i)/(s-p_B2_sD2(i));
        yf2_m_1 = yf2_m_1 + r_B2_sD2(i)/(s-p_B2_sD2(i));
    elseif r_B2_sD2(i)<0
        ye2_m_1 = ye2_m_1 - r_B2_sD2(i)/(s-p_B2_sD2(i));
        yf2_m_2 = yf2_m_2 - r_B2_sD2(i)/(s-p_B2_sD2(i));
    end
end

if size(k_B2_sD2, 1) ~= 0
    for i = 1:size(k_B2_sD2, 1)
        if k_B2_sD2(i)>0
            ye2_m_2 = ye2_m_2 + k_B2_sD2*s*(i-1);
            yf2_m_1 = yf2_m_1 + k_B2_sD2*s*(i-1);
        elseif k_B2_sD2(i)<0
            ye2_m_1 = ye2_m_1 - k_B2_sD2*s*(i-1);
            yf2_m_2 = yf2_m_2 - k_B2_sD2*s*(i-1);
        end
    end
end

if yf2_m_1 == 0
    yf2_m = vpa(simplify(s*ye2_m_2));
    ye2_m = vpa(simplify(s*yf2_m_2));
    ya2_m = vpa(simplify(s*ya2_m_2));
    yb2_m = vpa(simplify(s*yb2_m_2));
elseif yf2_m_2 == 0
    yf2_m = vpa(simplify(s*ye2_m_1));
    ye2_m = vpa(simplify(s*yf2_m_1));
    ya2_m = vpa(simplify(s*ya2_m_1));
    yb2_m = vpa(simplify(s*yb2_m_1));
else
    yf2_m = vpa(simplify(s*ye2_m_1));
    ye2_m = vpa(simplify(s*yf2_m_1));
    ya2_m = vpa(simplify(s*ya2_m_1));
    yb2_m = vpa(simplify(s*yb2_m_1));
end

yc2_m = yb2_m + yf2_m;
yd2_m = ya2_m + ye2_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% kuh F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degree_of_D_F1 = length(DEN1)-2; %1+1 for constant in DEN1
roots_of_DS_F1 = -1 * randperm((3*degree_of_D_F1),degree_of_D_F1);

DS_F1 = [1 -roots_of_DS_F1(1)];
for i = 2:degree_of_D_F1
    DS_temp_F1 = [1 -roots_of_DS_F1(i)];
    DS_F1 = conv(DS_F1, DS_temp_F1);
end

s = sym('s'); 
syms y11_F1 y12_F1 y21_F1 y22_F1
y11_F1 = 0; %also alpha=1
y21_F1 = -1 * (tf(NUM1,DS_F1));
y12_F1 = 0; %inintialize
y22_F1 = 0; %inintialize

NUM_bsd_F1 = DEN1;
DEN_bsd_F1 = conv([1 0],DS_F1);
bsd_F1 = tf(NUM_bsd_F1, DEN_bsd_F1);  %bsd is B/sD

[r_F1,p_F1,k_F1] = residue(NUM_bsd_F1,DEN_bsd_F1);
 for i=1:length(r_F1)
     if r_F1(i)>0
     y22_F1 = y22_F1 + tf(2*r_F1(i),[1 -p_F1(i)]);
     y12_F1 = y12_F1 + tf(abs(r_F1(i)),[1 -p_F1(i)]);
     elseif r_F1(i)<0
     y22_F1 = y22_F1 + tf(r_F1(i),[1 -p_F1(i)]);
     y12_F1 = y12_F1 + tf(2*abs(r_F1(i)),[1 -p_F1(i)]);
     end
 end
 if k_F1>0
     y22_F1 = y22_F1 + k_F1;
 elseif k_F1<0
     y12_F1 = y12_F1 + k_F1;
 end
 
[num_y12_F1,den_y12_F1] = tfdata(y12_F1);
num_y12_F1 = cell2mat(num_y12_F1);
den_y12_F1 = cell2mat(den_y12_F1);

[num_y22_F1,den_y22_F1] = tfdata(y22_F1);
num_y22_F1 = cell2mat(num_y22_F1);
den_y22_F1 = cell2mat(den_y22_F1);

new_den_y22_F1 =zeros(1,length(den_y22_F1)-1);
 for i=1:length(den_y22_F1)-1
     new_den_y22_F1(i) = den_y22_F1(i); %in new_den_y22_F1 we delete an s from den_y22_F1 to convert B/sD to B/D
 end
 y22_F1 = tf(num_y22_F1,new_den_y22_F1);
 
new_den_y12_F1 =zeros(1,length(den_y12_F1)-1);
 for i=1:length(den_y12_F1)-1
     new_den_y12_F1(i) = den_y12_F1(i);
 end
 y12_F1 = tf(num_y12_F1,new_den_y12_F1);

tf(num_y22_F1,new_den_y22_F1)
y22_F1_func_is_suitable_for_cauer1 = true;

%perform cauer 1 on y22_F1
%%if func has more than 20 k add to the dimension of kn_F1 in below line.
kn_F1 = zeros(1,20);
new_num_FLCs_F1 = num_y22_F1;
new_den_FLCs_F1 = new_den_y22_F1;
new_den_sFLCs_F1 = conv([1 0],new_den_FLCs_F1); %sFLCs stand for (1/s)*FLCs

sys_syms = poly2sym(new_num_FLCs_F1,s)/poly2sym(new_den_sFLCs_F1,s);
% limit sFLCs in s->inf
kn_F1(1) = limit(sys_syms, s, inf);

if kn_F1(1)<0
    disp("kuh method can not be performed with y22_F1")
    y22_F1_func_is_suitable_for_cauer1 = false;
end
num_FLCs_F1 = new_den_FLCs_F1;
den_FLCs_F1 = new_num_FLCs_F1 - (kn_F1(1)*new_den_sFLCs_F1);
den_sFLCs_F1 = conv([1 0], den_FLCs_F1);

%delete zeros at the begining of den_FLCs and den_sFLCs
counter1 = 0;
counter2 = 0;
for j=1:length(den_FLCs_F1)
    if den_FLCs_F1(j)==0
        counter1 = counter1 + 1;
    else
        break
    end
end

for jj=1:length(den_sFLCs_F1)
    if den_sFLCs_F1(jj)==0
        counter2 = counter2 + 1;
    else
        break
    end
end           
den_FLCs_F1 = den_FLCs_F1(counter1+1:end);
den_sFLCs_F1 = den_sFLCs_F1(counter2+1:end);


%while
i = 2;
   while (length(den_FLCs_F1)~=1) && y22_F1_func_is_suitable_for_cauer1
       if (length(num_FLCs_F1) == length(den_FLCs_F1)) && (rem(i,2)==0) %deconv process
           %disp('deconv started...')
           [q_F1,r_F1] = deconv(num_FLCs_F1,den_FLCs_F1);
           kn_F1(i) = q_F1;
           if kn_F1(i)<0
               disp("kuh method can not be performed with y22_F1")
               y22_F1_func_is_suitable_for_cauer1 = false;
           end
           num_FLCs_F1 = den_FLCs_F1;
           den_FLCs_F1 = r_F1;
           den_sFLCs_F1 = conv([1 0],den_FLCs_F1);
           
           % delete zeros in the begining of den_FLCs and den_sFLCs
           counter1 = 0;
           counter2 = 0;
           for j=1:length(den_FLCs_F1)
                if den_FLCs_F1(j)==0
                    counter1 = counter1 + 1;
                else
                    break
                end
            end

            for jj=1:length(den_sFLCs_F1)
                if den_sFLCs_F1(jj)==0
                    counter2 = counter2 + 1;
                else
                    break
                end
            end
            den_FLCs_F1 = den_FLCs_F1(counter1+1:end);
            den_sFLCs_F1 = den_sFLCs_F1(counter2+1:end);
                  
       elseif (length(num_FLCs_F1) == length(den_FLCs_F1)+1) && (rem(i,2)==1) %cauer1 process
                %disp("cauer1 preocess started...")
                new_num_FLCs_F1 = num_FLCs_F1;
                new_den_FLCs_F1 = den_FLCs_F1;
                new_den_sFLCs_F1 = den_sFLCs_F1;
                
                sys_syms = poly2sym(new_num_FLCs_F1,s)/poly2sym(new_den_sFLCs_F1,s);
                kn_F1(i) = limit(sys_syms, s, inf);
                
                if kn_F1(i)<0
                    disp("kuh method can not be performed with y22_F1")
                    y22_F1_func_is_suitable_for_cauer1 = false;
                end
                
                num_FLCs_F1 = new_den_FLCs_F1;
                den_FLCs_F1 = new_num_FLCs_F1 - (kn_F1(i).*new_den_sFLCs_F1);
                den_sFLCs_F1 = conv([1 0],den_FLCs_F1);
                
                %calculation Error in matlab
                num_FLCs_F1(abs(num_FLCs_F1)<0.000001) = 0; 
                den_FLCs_F1(abs(den_FLCs_F1)<0.000001) = 0;
                den_sFLCs_F1(abs(den_sFLCs_F1)<0.000001) = 0;
                
                %delete zeros at the begining of den_FLCs and den_sFLCs
                counter1 = 0;
                counter2 = 0;
                for j=1:length(den_FLCs_F1)
                    if den_FLCs_F1(j)==0
                        counter1 = counter1 + 1;
                    else
                        break
                    end
                end

                for jj=1:length(den_sFLCs_F1)
                    if den_sFLCs_F1(jj)==0
                        counter2 = counter2 + 1;
                    else
                        break
                    end
                end
           
                den_FLCs_F1 = den_FLCs_F1(counter1+1:end);
                den_sFLCs_F1 = den_sFLCs_F1(counter2+1:end);

       end
       i = i+1;
   end
   if y22_F1_func_is_suitable_for_cauer1
        kn_F1(i) = num_FLCs_F1(1)/den_FLCs_F1;
   end
    if kn_F1(i)<0
        disp("kuh method can not be performed with y22_F1")
        y22_F1_func_is_suitable_for_cauer1 = false;
    end
   if y22_F1_func_is_suitable_for_cauer1
        i = i+1;
        kn_F1(i) = num_FLCs_F1(2)/den_FLCs_F1;
   end
       if kn_F1(i)<0
            disp("kuh method can not be performed with y22_F1")
            y22_F1_func_is_suitable_for_cauer1 = false;
       end
       if ~y22_F1_func_is_suitable_for_cauer1
           kn_F1 = 0;
       end
       
   %delete unused zeros at the end of kn_F1
   for iii=1:length(kn_F1)
       if kn_F1(iii) == 0
           break
       end
   end
   kn_F1 = kn_F1(1:iii-1);
   
% seperation of kn_F1 to R and C
   R_y22_F1_kuh = zeros(1,floor(length(kn_F1)/2));
   C_F1_y22_F1_kuh = zeros(1,floor(length(kn_F1)/2));
   index_C = 1;
   index_R =1;
   for iter=1:length(kn_F1)
       if rem(iter,2) == 1
           C_F1_y22_F1_kuh(index_C) = kn_F1(iter);
           index_C = index_C + 1;
       elseif rem(iter,2) == 0
           R_y22_F1_kuh(index_R) = kn_F1(iter);
           index_R = index_R + 1;
       end
   end
   if y22_F1_func_is_suitable_for_cauer1
        kn_F1
        R_y22_F1_kuh
        C_F1_y22_F1_kuh
   end
   if ~y22_F1_func_is_suitable_for_cauer1
        disp("thus there is no answer for kuh method for F2")
   end
   Y_F1 = [y11_F1 y12_F1 ; y21_F1 y22_F1];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Kuh F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Active Design F2 (Kuh methode)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degree_of_D_F2 = length(DEN2)-2; %1+1 for constant in DEN2
roots_of_DS_F2 = -1 * randperm((3*degree_of_D_F2),degree_of_D_F2);

DS_F2 = [1 -roots_of_DS_F2(1)];
for i = 2:degree_of_D_F2
    DS_temp_F2 = [1 -roots_of_DS_F2(i)];
    DS_F2 = conv(DS_F2, DS_temp_F2);
end

s = sym('s'); 
syms y11_F2 y12_F2 y21_F2 y22_F2
y11_F2 = 0; %also alpha=1
y21_F2 = -1 * (tf(NUM2,DS_F1));
y12_F2 = 0; %inintialize
y22_F2 = 0; %inintialize

NUM_bsd_F2 = DEN2;
DEN_bsd_F2 = conv([1 0],DS_F2);
bsd_F2 = tf(NUM_bsd_F2, DEN_bsd_F2);  %bsd is B/sD

[r_F2,p_F2,k_F2] = residue(NUM_bsd_F2,DEN_bsd_F2);
 for i=1:length(r_F2)
     if r_F2(i)>0
     y22_F2 = y22_F2 + tf(2*r_F2(i),[1 -p_F2(i)]);
     y12_F2 = y12_F2 + tf(abs(r_F2(i)),[1 -p_F2(i)]);
     elseif r_F2(i)<0
     y22_F2 = y22_F2 + tf(r_F2(i),[1 -p_F2(i)]);
     y12_F2 = y12_F2 + tf(2*abs(r_F2(i)),[1 -p_F2(i)]);
     end
 end
 if k_F2>0
     y22_F2 = y22_F2 + k_F2;
 elseif k_F2<0
     y12_F2 = y12_F2 + k_F2;
 end
 
[num_y12_F2,den_y12_F2] = tfdata(y12_F2);
num_y12_F2 = cell2mat(num_y12_F2);
den_y12_F2 = cell2mat(den_y12_F2);

[num_y22_F2,den_y22_F2] = tfdata(y22_F2);
num_y22_F2 = cell2mat(num_y22_F2);
den_y22_F2 = cell2mat(den_y22_F2);

new_den_y22_F2 =zeros(1,length(den_y22_F2)-1);
 for i=1:length(den_y22_F2)-1
     new_den_y22_F2(i) = den_y22_F2(i);
 end
 y22_F2 = tf(num_y22_F2,new_den_y22_F2);
 
new_den_y12_F2 =zeros(1,length(den_y12_F2)-1);
 for i=1:length(den_y12_F2)-1
     new_den_y12_F2(i) = den_y12_F2(i);
 end
 y12_F2 = tf(num_y12_F2,new_den_y12_F2);
 
 
 
%perform cauer 1 on y22_F2
%if func has more than 20 k add to the dimension of kn_F2 in below line.
y22_F2_func_is_suitable_for_cauer1 = true;
kn_F2 = zeros(1,20);
new_num_FLCs_F2 = num_y22_F2;
new_den_FLCs_F2 = new_den_y22_F2;
new_den_sFLCs_F2 = conv([1 0],new_den_FLCs_F2); %sFLCs stand for (1/s)*FLCs

sys_syms = poly2sym(new_num_FLCs_F2,s)/poly2sym(new_den_sFLCs_F2,s);
% limit sFLCs in s->inf
kn_F2(1) = limit(sys_syms, s, inf);

if kn_F2(1)<0
    disp("kuh method can not be performed with y22_F2")
    y22_F2_func_is_suitable_for_cauer1 = false;
end
num_FLCs_F2 = new_den_FLCs_F2;
den_FLCs_F2 = new_num_FLCs_F2 - (kn_F2(1)*new_den_sFLCs_F2);
den_sFLCs_F2 = conv([1 0], den_FLCs_F2);

%delete zeros at the begining of den_FLCs and den_sFLCs
counter1 = 0;
counter2 = 0;
for j=1:length(den_FLCs_F2)
    if den_FLCs_F2(j)==0
        counter1 = counter1 + 1;
    else
        break
    end
end

for jj=1:length(den_sFLCs_F2)
    if den_sFLCs_F2(jj)==0
        counter2 = counter2 + 1;
    else
        break
    end
end           
den_FLCs_F2 = den_FLCs_F2(counter1+1:end);
den_sFLCs_F2 = den_sFLCs_F2(counter2+1:end);


%while
i = 2;
   while (length(den_FLCs_F2)~=1) && y22_F2_func_is_suitable_for_cauer1
       if (length(num_FLCs_F2) == length(den_FLCs_F2)) && (rem(i,2)==0) %deconv process
           %disp('deconv started...')
           [q_F2,r_F2] = deconv(num_FLCs_F2,den_FLCs_F2);
           kn_F2(i) = q_F2;
           if kn_F2(i)<0
               disp("kuh method can not be performed with y22_F2")
               y22_F2_func_is_suitable_for_cauer1 = false;
           end
           num_FLCs_F2 = den_FLCs_F2;
           den_FLCs_F2 = r_F2;
           den_sFLCs_F2 = conv([1 0],den_FLCs_F2);
           
           % delete zeros in the begining of den_FLCs and den_sFLCs
           counter1 = 0;
           counter2 = 0;
           for j=1:length(den_FLCs_F2)
                if den_FLCs_F2(j)==0
                    counter1 = counter1 + 1;
                else
                    break
                end
            end

            for jj=1:length(den_sFLCs_F2)
                if den_sFLCs_F2(jj)==0
                    counter2 = counter2 + 1;
                else
                    break
                end
            end
            den_FLCs_F2 = den_FLCs_F2(counter1+1:end);
            den_sFLCs_F2 = den_sFLCs_F2(counter2+1:end);
                  
       elseif (length(num_FLCs_F2) == length(den_FLCs_F2)+1) && (rem(i,2)==1) %cauer1 process
                %disp("cauer1 preocess started...")
                new_num_FLCs_F2 = num_FLCs_F2;
                new_den_FLCs_F2 = den_FLCs_F2;
                new_den_sFLCs_F2 = den_sFLCs_F2;
                
                sys_syms = poly2sym(new_num_FLCs_F2,s)/poly2sym(new_den_sFLCs_F2,s);
                kn_F2(i) = limit(sys_syms, s, inf);
                
                if kn_F2(i)<0
                    disp("kuh method can not be performed with y22_F2")
                    y22_F2_func_is_suitable_for_cauer1 = false;
                end
                
                num_FLCs_F2 = new_den_FLCs_F2;
                den_FLCs_F2 = new_num_FLCs_F2 - (kn_F2(i).*new_den_sFLCs_F2);
                den_sFLCs_F2 = conv([1 0],den_FLCs_F2);
                
                %calculation Error in matlab
                num_FLCs_F2(abs(num_FLCs_F2)<0.000001) = 0; 
                den_FLCs_F2(abs(den_FLCs_F2)<0.000001) = 0;
                den_sFLCs_F2(abs(den_sFLCs_F2)<0.000001) = 0;
                
                %delete zeros at the begining of den_FLCs and den_sFLCs
                counter1 = 0;
                counter2 = 0;
                for j=1:length(den_FLCs_F2)
                    if den_FLCs_F2(j)==0
                        counter1 = counter1 + 1;
                    else
                        break
                    end
                end

                for jj=1:length(den_sFLCs_F2)
                    if den_sFLCs_F2(jj)==0
                        counter2 = counter2 + 1;
                    else
                        break
                    end
                end
           
                den_FLCs_F2 = den_FLCs_F2(counter1+1:end);
                den_sFLCs_F2 = den_sFLCs_F2(counter2+1:end);

       end
       i = i+1;
   end
   if y22_F2_func_is_suitable_for_cauer1
        kn_F2(i) = num_FLCs_F2(1)/den_FLCs_F2;
   end
    if kn_F2(i)<0
        disp("kuh method can not be performed with y22_F2")
        y22_F2_func_is_suitable_for_cauer1 = false;
    end
   if y22_F2_func_is_suitable_for_cauer1
        i = i+1;
        kn_F2(i) = num_FLCs_F2(2)/den_FLCs_F2;
   end
       if kn_F2(i)<0
            disp("kuh method can not be performed with y22_F2")
            y22_F2_func_is_suitable_for_cauer1 = false;
       end
       if ~y22_F2_func_is_suitable_for_cauer1
           kn_F2 = 0;
       end
       
   %delete unused zeros at the end of kn_F2
   for iii=1:length(kn_F2)
       if kn_F2(iii) == 0
           break
       end
   end
   kn_F2 = kn_F2(1:iii-1);
   
   % seperation of kn_F2 to R and C
   R_y22_F2_kuh = zeros(1,floor(length(kn_F2)/2));
   C_F1_y22_F2_kuh = zeros(1,floor(length(kn_F2)/2));
   index_C = 1;
   index_R =1;
   for iter=1:length(kn_F2)
       if rem(iter,2) == 1
           C_F1_y22_F2_kuh(index_C) = kn_F2(iter);
           index_C = index_C + 1;
       elseif rem(iter,2) == 0
           R_y22_F2_kuh(index_R) = kn_F2(iter);
           index_R = index_R + 1;
       end
   end
   if y22_F2_func_is_suitable_for_cauer1
        kn_F2
        R_y22_F2_kuh
        C_F1_y22_F2_kuh
   end
   if ~y22_F2_func_is_suitable_for_cauer1
        disp("thus there is no answer for kuh method for F2")
   end
   Y_F2 = [y11_F2 y12_F2 ; y21_F2 y22_F2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% State Variable F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G1 = zeros(1, 2*N+1);
G1(1) = NUM1;
G1_final = G1;
for i = 2:2*N+1
    G1(i) = DEN1(i);
    G1_final(i) = Rg_denorm/G1(i);
end
R1_final = 1*Rg_denorm;
C1_final = 1/Rg_denorm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% State Variable F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G2 = zeros(1, 2*N+1);
G2(1) = NUM2;
G2_final = G2;
for i = 2:2*N+1
    G2(i) = DEN2(i);
    G2_final(i) = Rg_denorm/G2(i);
end
R2_final = 1*Rg_denorm;
C2_final = 1/Rg_denorm;