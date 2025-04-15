% Ikki, S., & Ahmed, M. H. (2007). 
% Performance analysis of cooperative diversity wireless networks over Nakagami-m fading channel. 
% IEEE Communications Letters, 11(4), 334–336.
% Fig.2

% Signal to Noise Ratio (SNR)
SNR1 = 0:1:25; 

N = 1e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent

ber_analytical1 = Lower_Bound_Func1(SNR1, A, N, 0.5, 1);
ber_analytical2 = Lower_Bound_Func1(SNR1, A, N, 0.5, 2);
ber_analytical3 = Lower_Bound_Func1(SNR1, A, N, 0.5, 3);
ber_simulation1 = Multipath_Simulation(SNR1, A, N, 0.5, 1);
ber_simulation2 = Multipath_Simulation(SNR1, A, N, 0.5, 2);
ber_simulation3 = Multipath_Simulation(SNR1, A, N, 0.5, 3);
ber_analytical4 = Lower_Bound_Func2(SNR1, A, N, 0.5, 1);
ber_analytical5 = Lower_Bound_Func2(SNR1, A, N, 0.5, 2);
ber_analytical6 = Lower_Bound_Func2(SNR1, A, N, 0.5, 3);


semilogy(SNR1, ber_analytical1, '--o red', SNR1, ber_simulation1, '-o red', SNR1, ber_analytical2, '--diamond blue',  SNR1, ber_simulation2, '-diamond blue', SNR1, ber_analytical3, '--^ green', SNR1, ber_simulation3, '-^ green', SNR1, ber_analytical4, ':o red', SNR1, ber_analytical5, ':diamond blue', SNR1, ber_analytical6, ':^ green ');
ylim([10^-7,10^0])
xlim([0,25])
title('BER: BPSK, Results of Several Number of Cooperative Paths under Nakagami-m Fading')
legend('Lower Bound M = 1 (Analytical eq.9)', 'Exact M=1 (Simulation)','Lower Bound M = 2 (Analytical eq.9)', 'Exact M=2 (Simulation)','Lower Bound M = 3 (Analytical eq.9)', 'Exact M=3 (Simulation)', 'Lower Bound M = 1 (Analytical eq.[4])', 'Lower Bound M = 2 (Analytical eq.[4])', 'Lower Bound M = 3 (Analytical eq.[4])');
xlabel('Es/N0 [dB]');
ylabel('Bit Error Rate (BER)');
grid on


function [ber_simulation] = Multipath_Simulation(SNR, A, N, m, M_path)
    absolute_amplitude = A; % Maximum amplitude of the signal
    M = 2; % BPSK 

    % Energy definitions
    Eb = (absolute_amplitude ^ 2) / 2;
    Es = Eb;

    % Voltage of the sent signals
    voltages = zeros(1,M);
    for i = 1:M
        theta = (i-1) *  180;
      
        voltages(i) = -Es * cosd(theta);
    end

    % Generate random symbols
    s_n = randsrc(1, N, voltages);

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Simulated error rate arrays
    ber_simulation = zeros(1, length(SNR));

    % Nakagami-m Distribution Characteristics
    h_sr_dist = makedist('Nakagami', 'mu', m, 'omega', 0.75);
    h_sd_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
    h_rd_dist = makedist('Nakagami', 'mu', m, 'omega', 0.5);
    
    for i = 1:length(SNR)
        z_n_cont = zeros(M_path, N);
        for j = 1:M_path
            h_sr = sqrt(1/2) .* (random(h_sr_dist, 1, N) + 1i  * random(h_sr_dist, 1, N));        
            h_rd = sqrt(1/2) .* (random(h_rd_dist, 1, N) + 1i  * random(h_rd_dist, 1, N));

            % Power Calculations
            energy_of_received_symbol_sr = mean(abs(h_sr).^2) * Es;
            N0_sr = energy_of_received_symbol_sr / (linear_SNR(i));
            energy_of_received_symbol_rd = mean(abs(h_rd).^2) * Es;
            N0_rd = energy_of_received_symbol_rd / (linear_SNR(i));
        
            % Phase 1 of Indirect Link
            x_r = h_sr .* s_n + sqrt(N0_sr / 2) .* (randn(1,N) + 1i  * randn(1,N));

            % Phase 2 of Indirect Link
            G = sqrt((N0_sr + abs(h_sr).^2));
            s_r = x_r ./ G;
            y_d = h_rd .* s_r + sqrt(N0_rd / 2) .* (randn(1,N) + 1i  * randn(1,N));

            % Linked Channel Properties
            linked_channel_h =  h_sr .* h_rd .* G;
            
            linked_channel_variance = (abs(h_rd).^2 .* N0_rd .* (N0_sr + abs(h_sr).^2).^(-1)) + N0_rd;
                        
            z_n_cont(j,:) = conj(linked_channel_h ./ linked_channel_variance) .* y_d ./ linked_channel_variance;
        end
               
        h_sd = sqrt(1/2) .* (random(h_sd_dist, 1, N) + 1i * random(h_sd_dist, 1, N));

        energy_of_received_symbol_sd = mean(abs(h_sd).^2) * Es;
        N0_sd = energy_of_received_symbol_sd / (linear_SNR(i));
        
        % Direct Link
        x_d = h_sd .* s_n + sqrt(N0_sd / 2) .* (randn(1,N) + 1i  * randn(1,N));
        
        % Decision Array
        decision = zeros(1, N);

        % Decision Variable

        if (M_path == 1)
            z_n = (conj(h_sd) .* x_d)  + z_n_cont;
        else
            z_n = (conj(h_sd) .* x_d) + sum(z_n_cont, 1);
        end
       
        % Decision Circuit
        decision(z_n > 0) = voltages(2);
        decision(z_n < 0) = voltages(1);

        % SER Simulation
        ber_simulation(i) = sum(decision ~= s_n) / N;
        disp(i);
    
    end
end

function [ber_analytical] = Lower_Bound_Func1(SNR, A, N, m, M_path)
    absolute_amplitude = A; % Maximum amplitude of the signal 

    % Energy definitions
    Eb = (absolute_amplitude ^ 2) / 2;
    Es = Eb;

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Analytical BER array
    ber_analytical = zeros(1, length(SNR));

    for i = 1:length(SNR)
        % Nakagami-m Fading Channels
        h_sd_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
        h_sr_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
        h_rd_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
        
        h_sr = sqrt(0.75/2) .* (random(h_sr_dist,1, N) + 1i * random(h_sr_dist,1, N));
        h_sd = sqrt(1/2) .* (random(h_sd_dist,1, N) + 1i * random(h_sd_dist,1, N));
        h_rd = sqrt(0.5/2) .* (random(h_rd_dist,1, N) + 1i * random(h_rd_dist,1, N));

        % Integral range and dx
        theta_range = linspace(pi/90, pi/2, 10);
        dx = (pi/2 - pi/90) / (length(theta_range) - 1);

        % Energy, Power & SNR Calculations
        energy_of_received_symbol_sr = mean(abs(h_sr).^2) * Es;
        N0_sr = energy_of_received_symbol_sr / (linear_SNR(i));
        avg_ins_snr_sr = mean(abs(h_sr).^2) * Es / N0_sr;
        
        energy_of_received_symbol_sd = mean(abs(h_sd).^2) * Es;
        N0_sd = energy_of_received_symbol_sd / (linear_SNR(i));
        avg_ins_snr_sd = mean(abs(h_sd).^2) * Es / N0_sd;
        mgf_sd = (1 + avg_ins_snr_sd / m .* sin(theta_range).^(-2) ).^(-m);

        energy_of_received_symbol_rd = mean(abs(h_rd).^2) * Es;
        N0_rd = energy_of_received_symbol_rd / (linear_SNR(i));
        avg_ins_snr_rd = mean(abs(h_rd).^2) * Es / N0_rd;

        % MGF Functions for S-->R & R-->D Links
        func_part1 = (m / avg_ins_snr_sr)^m * (m / avg_ins_snr_rd)^m * (gamma(2*m) / (gamma(m)^2));
        func_part2 = (1/m) .* hypergeom([1 2*m], (m+1), (m / avg_ins_snr_sr + sin(theta_range).^(-2)) ./ (m / avg_ins_snr_sr + m / avg_ins_snr_rd + sin(theta_range).^(-2)));
        func_part3 = (1/m) .* hypergeom([1 2*m], (m+1), (m / avg_ins_snr_rd + sin(theta_range).^(-2)) ./ (m / avg_ins_snr_sr + m / avg_ins_snr_rd + sin(theta_range).^(-2)));
        func_part4 = (1 ./ ((m / avg_ins_snr_rd) + (m / avg_ins_snr_sr) + sin(theta_range).^(-2)).^(2*m));
               
        mgf_func = func_part1 .* (func_part2 + func_part3) .* func_part4;
        % Lower Bound MGF Function for the system
        mgf_bottom = mgf_sd .* (mgf_func.^M_path);

        %BER Simulation
        ber_analytical(i) = (1/pi) .* sum(mgf_bottom) * dx;
        disp(i);
    end
end

function [ber_analytical] = Lower_Bound_Func2(SNR, A, N, m, M_path)
    absolute_amplitude = A; % Maximum amplitude of the signal
    M = 2; % BPSK 

    % Energy definitions
    Eb = (absolute_amplitude ^ 2) / 2;
    Es = Eb;

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Analytical BER array
    ber_analytical = zeros(1, length(SNR));

    % Distribution of fading channels
    h_sd_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
    h_sr_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
    h_rd_dist = makedist('Nakagami', 'mu', m, 'omega', 1);
    
    for i = 1:length(SNR)
        h_sr = sqrt(0.75/2) .* (random(h_sr_dist,1, N) + 1i * random(h_sr_dist,1, N));
        h_sd = sqrt(1/2) .* (random(h_sd_dist,1, N) + 1i * random(h_sd_dist,1, N));
        h_rd = sqrt(0.5/2) .* (random(h_rd_dist,1, N) + 1i * random(h_rd_dist,1, N));
    
        energy_of_received_symbol_sr = mean(abs(h_sr).^2) * Es;
        N0_sr = energy_of_received_symbol_sr / (linear_SNR(i));
        
        energy_of_received_symbol_sd = mean(abs(h_sd).^2) * Es;
        N0_sd = energy_of_received_symbol_sd / (linear_SNR(i));

        energy_of_received_symbol_rd = mean(abs(h_rd).^2) * Es;
        N0_rd = energy_of_received_symbol_rd / (linear_SNR(i));

        x = linspace(0.6, 10, 1000);
        p_gamma_sd = gampdf(x, m, linear_SNR(i) / m);        
        p_gamma_rd = gampdf(x, m, linear_SNR(i) / m);
        p_gamma_sr = gampdf(x, m, linear_SNR(i) / m);
        disp(trapz(x, p_gamma_rd))
        addition_term = 1;

        for j = 1:M_path+1
            addition_term = addition_term * (p_gamma_sr(1) + p_gamma_rd(1));   
        end

        % Calculation of C(M)
        numerative = 1;
        denumerative = ((2* factorial(M_path+1) * 1 ^(M_path+1)));
        for k = 1 : M_path+1
            numerative = numerative * (2*k - 1);  
        end
        C_M = numerative / denumerative;

        
        Pe = C_M * p_gamma_sd(1) * addition_term;
        ber_analytical(i) = Pe; 
    
    end
    




end

