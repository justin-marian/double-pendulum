# Double Pendulum

## Description

The Double Pendulum is a `MATLAB` code that simulates the **motion of three double pendulum systems** with different physical parameters. It calculates the motion of each pendulum using recurrence relations and then visualizes the motion using graphical simulation.

The dynamics of the double pendulum are described using **Lagrangian mechanics** and **Runge-Kutta (RK4)**, which provides a more elegant and efficient way to derive the equations of motion for complex systems like the double pendulum.

<div style="display: flex; justify-content: space-between;">
    <img src="/images/three_double_pendulums.png" alt="VELOCITY_POSITION_CURVE" style="width: 49%; height: auto;" />
    <img src="/images/runge_kutta_4ode.png" alt="RUNGE_KUTTA_4ODE" style="width: 49%; height: 400px;" />
</div>

## Physical Parameters

- **Pendulum Parameters:** Each pendulum system is defined by its lengths (`L1` and `L2`) and masses (`m1` and `m2`). These parameters determine the geometry and inertia of the pendulum systems, influencing their motion.
- **Initial Conditions:** Initial angles (`theta1` and `theta2`) and initial angular velocities (`OM1` and `OM2`)  specify the starting configuration and motion of each pendulum system. These conditions play a crucial role in determining the subsequent motion and behavior of the systems.

## Simulation Settings

- **Time Settings:** The simulation duration is determined by the characteristic time of the double pendulum motion. The simulation time is discretized into a specified number of time steps.
- **Recurrence Cycle:** Recurrence relations are used to calculate the motion of each pendulum system.
- **Energy Calculation:** The kinetic and potential energies of each pendulum system are calculated at each time simulation step (`1cs`).
- **Graphical Simulation:** The motion of the pendulums is visualized using a graphical simulation.
- **Tracebacks:** The pendulum positions and bottom objects' positions are also plotted during the simulation.

## Lagrangian Equations

The **Lagrangian** ($\mathcal{L}$) of the double pendulum is given by the difference between the **kinetic energy** (**$T$**) and the **potential energy** (**$U$**) of the system:

1. **Kinetic Energy (T):**
   \[ T = \frac{1}{2} m_1 L_1^2 \dot{\theta}_1^2 + \frac{1}{2} m_2 (L_1^2 \dot{\theta}_1^2 + L_2^2 \dot{\theta}_2^2 + 2 L_1 L_2 \dot{\theta}_1 \dot{\theta}_2 \cos(\theta_1 - \theta_2)) \]

2. **Potential Energy (U):**
   \[ U = -m_1 g L_1 \cos(\theta_1) - m_2 g (L_1 \cos(\theta_1) + L_2 \cos(\theta_2)) \]

The total energy of the double pendulum system, given by the ***sum of kinetic energy (T) and potential energy (U), must remain constant over time***. In theory, the total energy should be conserved throughout the simulation due to the **ideal** nature of the system, which assumes the *absence of friction or aerodynamic resistance*.

### Equations of Motion

The equations of motion for the double pendulum are derived by applying the **Euler-Lagrange** equations to the Lagrangian. These equations are ***second-order ordinary differential equations (ODE2)*** that describe the evolution of the angles $\theta_1$ and $\theta_2$ over time.

The **Euler-Lagrange** equations for the double pendulum are given by:

\[
\frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{\theta}_1}\right) - \frac{\partial \mathcal{L}}{\partial \theta_1} = 0
\]

\[
\frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{\theta}_2}\right) - \frac{\partial \mathcal{L}}{\partial \theta_2} = 0
\]

## Runge-Kutta (ODE4)

The **Runge-Kutta** is used for solving ordinary differential equations. Specifically, the ***fourth-order variant, RK4***, is widely used. It estimates the solution at the next time step based on the solution and its derivative at the current time step. By iteratively applying this process, RK4 generates a numerical solution to the ODE.

- **RK4 Process**: Estimates the solution at the next time step based on the current solution and its derivative. This iterative process generates a numerical solution to the ODE.
- **Discretization of Time Domain**: RK4 discretizes the time domain into small intervals, enabling the computation of angular positions and velocities of the pendulum at each step.

### Algorithm Runge-Kutta

This pseudocode outlines the implementation of the **Runge-Kutta (RK4)** method for numerically solving ordinary differential equations:

- **$k_1 \gets h \cdot f(t_n, y_n)$:** Computes the slope at $(t_n, y_n)$, scaled by the step size $h$.
- **$k_2 \gets h \cdot f(t_n + \frac{h}{2}, y_n + \frac{k_1}{2})$:** Computes the slope at the midpoint, using $k_1$ adjusted state.
- **$k_3 \gets h \cdot f(t_n + \frac{h}{2}, y_n + \frac{k_2}{2})$:** Computes slope at midpoint, using $k_2$ adjusted state.
- **$k_4 \gets h \cdot f(t_n + h, y_n + k_3)$:** Computes slope at the next step, using $k_3$ adjusted state.
- **$y_{n+1} \gets y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$:** Updates state using weighted average of slopes.
- **$t_{n+1} \gets t_n + h$:** Updates time to the next step.

## Numerical aproximations

Due to numerical approximations or errors in the simulation, there might be a small deviation in the total energy, typically within the range of **0.01 to 0.1**. Monitoring this deviation can help ensure the accuracy and stability of the simulation.
