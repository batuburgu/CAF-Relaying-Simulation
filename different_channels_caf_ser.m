% Signal to Noise Ratio (SNR)
SNR1 = 0:1:30; 
N = 1e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent

ber_simulation0 = Multipath_Simulation(SNR1, A, N, 1, 1, 1.6, 0);
ber_simulation1 = Multipath_Simulation(SNR1, A, N, 1, 1, 1.6, 1);
ber_simulation2 = Multipath_Simulation(SNR1, A, N, 1, 1, 1.6, 2);
ber_simulation3 = Multipath_Simulation(SNR1, A, N, 1, 1, 1.6, 3);
ber_simulation4 = Multipath_Simulation(SNR1, A, N, 1, 1, 1.6, 4);
ber_analytical0 = Theoretical_Result(SNR1, 1, 1, 1.6, 0);
ber_analytical1 = Theoretical_Result(SNR1, 1, 1, 1.6, 1);
ber_analytical2 = Theoretical_Result(SNR1, 1, 1, 1.6, 2);
ber_analytical3 = Theoretical_Result(SNR1, 1, 1, 1.6, 3);
ber_analytical4 = Theoretical_Result(SNR1, 1, 1, 1.6, 4);

semilogy(SNR1, ber_simulation0, '-square', SNR1, ber_simulation1, '-o', SNR1, ber_simulation2, '-diamond', SNR1, ber_simulation3, '-^', SNR1, ber_simulation4, '-x', SNR1, ber_analytical0, '--square', SNR1, ber_analytical1, '--o', SNR1, ber_analytical2, '--diamond', SNR1, ber_analytical3, '--^', SNR1, ber_analytical4, '--x')
ylim([10^-5,10^0])
xlim([0,30])
title('BER: BPSK, Results of Several Number of Cooperative Paths under Nakagami-m Fading Channels with Different m Parameters')
legend( 'Simulation (M = 0)', 'Simulation (M = 1)', 'Simulation (M = 2)', 'Simulation (M = 3)', 'Simulation (M = 4)', 'Analytical (M = 0)', 'Analytical (M = 1)', 'Analytical (M = 2)', 'Analytical (M = 3)', 'Analytical (M = 4)');
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
    linear_SNR = 10.^(SNR / 10) ./ (M_path + 1);

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
            y_d = h_rd .* s_r + sqrt(N0_rd / 2) .* (randn(1,N) + 1i  * randn(1,N));

            % Linked Channel Properties
            linked_channel_h =  h_sr .* h_rd ./ G;         
            linked_channel_variance = (abs(h_rd ./ G).^2 + 1) * N0_rd;
                        
            z_n_cont(j,:) = conj(linked_channel_h ./ sqrt(linked_channel_variance)) .* y_d ./ sqrt(linked_channel_variance); 
        end
        % Direct Channel Communication   
        h_sd = random(h_sd_dist, 1, N);
        energy_of_received_symbol_sd = mean(abs(h_sd).^2) * Es;
        N0_sd = energy_of_received_symbol_sd ./ linear_SNR(i);

        x_d = h_sd .* s_n + sqrt(N0_sd / 2) .* (randn(1,N) + 1i  * randn(1,N));
        
        % Decision Array
        decision = zeros(1, N);

        % Decision Variable
        if (M_path == 0)
            z_n = x_d;
        elseif (M_path == 1)
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

function [ber_analytical] = Theoretical_Result(SNR, m_sd, m_sr, m_rd, M_path)
    linear_SNR = 10.^(SNR / 10) / (M_path + 1);
    ber_analytical = zeros(1, length(SNR));
    theta_range = linspace(1e-6, pi/2, 10);

    for i = 1:length(SNR) 
        % Indirect Channel
        alpha = (m_sr/linear_SNR(i)) + (m_rd/linear_SNR(i)) + (1./(sin(theta_range).^2));
        beta = 2 * sqrt(m_sr * m_rd / (linear_SNR(i)^2) );
        mu = m_sr + m_rd;
        
        Il_part1 = (2 / (gamma(m_sr) * gamma(m_rd))) * ((m_sr/linear_SNR(i))^m_sr) * ((m_rd/linear_SNR(i))^m_rd);

        summation_part = 0;
        for k = 0:floor(m_sr + 2*m_rd)- 1
            v = abs(m_rd - k);
            
            Il_part2 = (m_sr * linear_SNR(i) / (m_rd * linear_SNR(i))).^((m_rd - k)/2); 
            Il_part3 = sqrt(pi) * (2*beta)^v .* gamma(mu + v) * gamma(mu - v) ./ (((alpha + beta).^(mu + v)) .* gamma(mu + 1/2));
            Il_part4 = (hypergeom([(mu + v) (v + 1/2)], (mu + 1/2), ((alpha - beta) ./ (alpha + beta))));

            summation_part = summation_part + A_func(m_sr + m_rd, k) .* Il_part2 .* Il_part3 .* Il_part4;

        end
        Il = (Il_part1 .* sum(summation_part,1));  
         
        % Direct Channel 
        I0 = ((m_sd .* sin(theta_range).^2) ./ (linear_SNR(i) + m_sd .* sin(theta_range).^2)).^m_sd;

        ber_analytical(i) = (1/pi) .* trapz(theta_range, I0 .* (Il.^M_path));

    end   
end

function [result] = A_func(q,k)
    if k == 0
        result = 1;
    else
        result = 1;
        for i = 0:(k-1)
            result = result * (q - i);
        end
        result = result / factorial(k);
    end
end
