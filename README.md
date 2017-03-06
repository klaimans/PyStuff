# PyStuff

# The Harmonic Interaction Model

### The N-body probabilty distribution function (PDF)

\begin{align} 
P_m &\equiv P(r_1,r_2,\ldots,r_m) = \left|\Psi\left(r_1,r_2,\ldots,r_N\right)\right|^2 \\
&= \left(\frac{\delta_N}{\pi}\right)^{\frac{dm}{2}}\left(\frac{N\omega}{(N-m)\omega +m\delta_N}\right)^{\frac{d}{2}}\exp\left(-\frac{1}{N}\left(\omega+(N-1)\delta_N\right)\sum_{i=1}^mr_i^2-\frac{2}{N}
\sum_{i<j}^mr_ir_j-C_m\left(\sum_{i=1}^mr_i\right)^2\right)
\end{align}


#### where
\begin{equation*}
C_m = -\frac{1}{N}\frac{(N-m)(\omega-\delta_N)^2}{(N-m)\omega + m \delta_N} 
\end{equation*}

### Probability "Chain Rule"

\begin{equation*}
P(r_1,r_2,\ldots,r_m) = P(r_1)P(r_2|r_1)\cdots P(r_m|r_1,r_2,\ldots,r_{m-1})
\end{equation*}

### The conditional probability

\begin{align} 
P(r_m|r_1,r_2,\ldots,r_{m-1}) &= \frac{P(r_1,r_2,\ldots,r_m)}{P(r_1,r_2,\ldots,r_{m-1})}= \frac{P_m}{P_{m-1}}\\ &= \left(\frac{\delta_N}{\pi}\right)^{\frac{d}{2}}\left(1 + \frac{\omega-\delta_N}{(N-m)\omega +m\delta_N}\right)^{\frac{d}{2}}\exp\left(-\frac{1}{N}(\omega+(N-1)\delta_N)r_m^2 -\frac{1}{N}r_m\sum_{i=1}^{m-1}r_i\right)\\
&\times\exp\left( -C_mr_m^2-2C_mr_m\sum_{i=1}^{m-1}r_i-\left(\sum_{i=1}^{m-1}r_i\right)^2\left(C_m-C_{m-1}\right) \right) \\
&= \frac{1}{\sqrt{2\pi\sigma_m^2}}\exp\left(-\frac{\left(r_m-\mu_m\right)^2}{2\sigma_m^2} \right)
\end{align}


#### where
\begin{equation*}
\mu_m = \frac{\omega-\delta_N}{(m-1)\delta_N+(1-m+N)\omega}\sum_{i=1}^{m-1}r_i 
\end{equation*}
\begin{equation*}
\sigma_m^2 = \frac{m(\delta_N-\omega)+N\omega}{2\delta_N\left((m-1)\delta_N+(1-m+N)\omega\right)} 
\end{equation*}

### The Boson Density

\begin{align} 
\rho_{1}(r,r^\prime) &= N \left(\frac{\delta_NN\omega/\pi}{(N-1)\omega+\delta_N}\right)^{d/2}\exp\left(-a_1\left(r^2 + r^{\prime 2} \right)+a_2rr^\prime\right) \\
&= N \left(\frac{\delta_NN\omega/\pi}{(N-1)\omega+\delta_N}\right)^{d/2}\exp\left(-\frac{1}{4}\left(a^2+b^2\right)\left(r^2 + r^{\prime 2} \right)-\frac{1}{2}\left(a^2-b^2\right)rr^\prime\right)
\end{align}


