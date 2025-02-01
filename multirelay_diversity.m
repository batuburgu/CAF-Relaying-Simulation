% Signal to Noise Ratio (SNR)
SNR1 = 0:1:24; 
N = 5e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent

ber_simulation1 = Multipath_Simulation(SNR1, A, N, 1, 1, 1, 1);
ber_simulation2 = Multipath_Simulation(SNR1, A, N, 1, 1, 1, 2);
ber_analytical0 = Theoretical_Result(SNR1, 0);
ber_analytical1 = Theoretical_Result(SNR1, 1);
ber_analytical2 = Theoretical_Result(SNR1, 2);
ber_analytical3 = Theoretical_Result(SNR1, 3);


semilogy(SNR1, ber_simulation1, '--', SNR1, ber_simulation2, '--diamond', SNR1, ber_analytical0, '-o', SNR1, ber_analytical1, SNR1, ber_analytical2, '-diamond', SNR1, ber_analytical3, '-square');
ylim([10^-8,10^0])
xlim([0,24])
title('BER: BPSK, Diversity Advantage of Multirelay Network')
legend('Simulation (K = 1)', 'Simulation (K = 2)', 'Analytical (K = 0)', 'Analytical (K = 1)', 'Analytical (K = 2)', 'Analytical (K = 3)');
xlabel('Es/N0 [dB]');
ylabel('Bit Error Rate (BER)');
grid on

function [ber_simulation] = Multipath_Simulation(SNR, A, N, m_sd, m_sr, m_rd, M_path)
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
    linear_SNR = 10.^(SNR / 10) / (M_path + 1);
    
    % Simulated error rate arrays
    ber_simulation = zeros(1, length(SNR));

    % Nakagami-m Distribution Characteristics
    h_sr_dist = makedist('Nakagami', 'mu', m_sr, 'omega', 1);
    h_sd_dist = makedist('Nakagami', 'mu', m_sd, 'omega', 1);
    h_rd_dist = makedist('Nakagami', 'mu', m_rd, 'omega', 1);
    
    for i = 1:length(SNR)
        z_n_cont = zeros(M_path, N);
        for j = 1:M_path
            h_sr = random(h_sr_dist, 1, N);        
            h_rd = random(h_rd_dist, 1, N);

            % Power Calculations
            energy_of_received_symbol_sr = mean(abs(h_sr).^2) * Es;
            N0_sr = energy_of_received_symbol_sr ./ linear_SNR(i);
            energy_of_received_symbol_rd = mean(abs(h_rd).^2) * Es;
            N0_rd = energy_of_received_symbol_rd ./ linear_SNR(i);
        
            % Phase 1
            x_r = h_sr .* s_n + sqrt(N0_sr / 2) .* (randn(1, N) + 1i * randn(1, N));

            % Phase 2
            G = sqrt((N0_sr + abs(h_sr).^2));
            s_r = x_r ./ G;
            y_d = h_rd .*  s_r + sqrt(N0_rd / 2) .* (randn(1,N) + 1i  * randn(1,N));

            % Linked Channel Properties
            linked_channel_h =  h_sr .* h_rd ./ G;         
            linked_channel_variance = (abs(h_rd ./ G).^2 + 1) * N0_rd;
                        
            z_n_cont(j,:) = conj(linked_channel_h ./ sqrt(linked_channel_variance)) .* y_d ./ sqrt(linked_channel_variance);  
        end
        % Direct Channel Communication   
        h_sd = random(h_sd_dist, 1, N);
        energy_of_received_symbol_sd = mean(abs(h_sd).^2) * Es;
        N0_sd = energy_of_received_symbol_sd ./ linear_SNR(i);

        x_d = h_sd .* s_n + sqrt(N0_sd  / 2) .* (randn(1,N) + 1i  * randn(1,N));
        
        % Decision Array
        decision = zeros(1, N);

        % Decision Variable
        if (M_path == 1)
            z_n = (conj(h_sd / sqrt(N0_sd)) .* x_d / sqrt(N0_sd)) + z_n_cont;
        else
            z_n = (conj(h_sd / sqrt(N0_sd)) .* x_d / sqrt(N0_sd)) + sum(z_n_cont, 1);
        end
       
        % Decision Circuit
        decision(z_n > 0) = voltages(2);
        decision(z_n < 0) = voltages(1);

        % SER Simulation
        ber_simulation(i) = sum(decision ~= s_n) / N;
        disp(i);
    
    end
end

function [ber_analytical] = Theoretical_Result(SNR, M_path)
    % Linear SNR
    linear_SNR = 10.^(SNR / 10) ./ (M_path + 1);
    
    % Results Array
    ber_analytical = zeros(1, length(SNR));

    % Analytical Calculations
    if M_path == 0
        for i = 1:length(SNR)
            Ps_part1 = @(theta) (1 + linear_SNR(i) ./ (sin(theta).^2)).^(-1);
            ber_analytical(i) = (1/pi) .* integral(Ps_part1, 0, pi/2);
        end

    else
        for i=1:length(SNR)
            gamma_sigma_k = 2 * linear_SNR(i);
            gamma_p_k = (linear_SNR(i))^2;
            Ik = @(theta) gamma_sigma_k + (gamma_p_k ./ (sin(theta).^2));

            Ps_part1 = @(theta) (1 + (linear_SNR(i) ./ (sin(theta).^2))).^(-1);

            Ps_part211 = @(theta) (4 .* gamma_p_k .* (Ik(theta) - gamma_sigma_k) ./ sqrt(Ik(theta).^2 - 4 * gamma_p_k));
            Ps_part212 = @(theta) log((Ik(theta) + sqrt(Ik(theta).^2 - 4 * gamma_p_k)) ./ (2 * sqrt(gamma_p_k)));
            Ps_part213 = @(theta) (gamma_sigma_k .* Ik(theta)) - (4 * gamma_p_k);
            Ps_part21 = @(theta) (Ps_part211(theta) .* Ps_part212(theta) + Ps_part213(theta));

            Ps_part22 = @(theta) (Ik(theta).^2 - 4 * gamma_p_k);
            
            Ps_part2 = @(theta) (Ps_part21(theta) ./ Ps_part22(theta));

            Ps = @(theta) Ps_part1(theta) .* Ps_part2(theta).^(M_path);

            ber_analytical(i) = (1/pi) .* integral(Ps, 0, pi/2);
        end
    end 


end
