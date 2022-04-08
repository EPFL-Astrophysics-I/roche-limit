// An abstract class for integrating a system of ODEs
public abstract class Integrator
{
    private int numEquations;
    private double[] store;
    private double[] k1;
    private double[] k2;
    private double[] k3;
    private double[] k4;
    private double[] ym1;
    private double[] ym2;
    private double[] ym3;
    private double[] P;
    private double[] dm1;
    private double[] dm2;
    private double[] dm3;
    private double[] dp1;
    private int abmSteps = 0;
    private double abmRms2;

    public Integrator()
    {
        Init(1);
    }

    // Allocate memory for all storage arrays and set number of equations
    public void Init(int numEquations)
    {
        // set up temp arrays
        this.numEquations = numEquations;
        store = new double[numEquations];
        k1 = new double[numEquations];
        k2 = new double[numEquations];
        k3 = new double[numEquations];
        k4 = new double[numEquations];
        ym1 = new double[numEquations];
        ym2 = new double[numEquations];
        ym3 = new double[numEquations];
        P = new double[numEquations];
        dm1 = new double[numEquations];
        dm2 = new double[numEquations];
        dm3 = new double[numEquations];
        dp1 = new double[numEquations];
        abmSteps = 0;
    }

    /// <summary>
    /// Abstract void, override this method to set the ODEs to be
    /// integrated.
    /// </summary>
    /// <param name="x">The values being integrated.</param>
    /// <param name="xdot">The derivatives being calculated.</param>
    public abstract void RatesOfChange(double[] x, double[] xdot, double t);

    /// <summary>
    /// Step forward using Euler's method
    /// </summary>
    /// <param name="x">The values being integrated.</param>
    /// <param name="h">The time step.</param>
    public void EulerStep(double[] x, double t, double h)
    {
        RatesOfChange(x, k1, t);
        for (int i = 0; i < numEquations; i++)
        {
            x[i] += k1[i] * h;
        }
    }

    /// <summary>
    /// Step forward using 4th order Runge Kutta method
    /// </summary>
    /// <param name="x">The values being integrated.</param>
    /// <param name="h">The time step.</param>
    public double RK4Step(double[] x, double t, double h)
    {
        RatesOfChange(x, k1, t);
        for (int i = 0; i < numEquations; i++)
        {
            store[i] = x[i] + k1[i] * h / 2.0;
        }
        RatesOfChange(store, k2, t);
        for (int i = 0; i < numEquations; i++)
        {
            store[i] = x[i] + k2[i] * h / 2.0;
        }
        RatesOfChange(store, k3, t);
        for (int i = 0; i < numEquations; i++)
        {
            store[i] = x[i] + k3[i] * h;
        }
        RatesOfChange(store, k4, t);
        for (int i = 0; i < numEquations; i++)
        {
            x[i] = x[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * h / 6.0;
        }
        return t + h;
    }


    /**
	 * Calculates a single step using Adams Bashforth Moulton,
	 * 
	 * @param x Array of values being integrated.
	 * @param t Time at which step begins
	 * @param h Duration of step
	 * @return Error prediction at end of step
	 */
    public double abmStep(double[] x, double t, double h)
    {
        abmRms2 = 0.0;
        if (abmSteps == 0)
        {
            for (int i = 0; i < x.Length; i++)
            {
                ym3[i] = x[i];
                ym2[i] = x[i];
            }
            RatesOfChange(dm3, ym3, t);
            t = RK4Step(ym2, t, h);
            RatesOfChange(dm2, ym2, t);
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = ym2[i];
            }
            abmSteps += 1;
            return 1.0;
        }
        else if (abmSteps == 1)
        {
            for (int i = 0; i < x.Length; i++)
            {
                ym1[i] = ym2[i];
            }
            t = RK4Step(ym1, t, h);
            RatesOfChange(dm1, ym1, t);
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = ym1[i];
            }
            abmSteps += 1;
            return 1.0;
        }
        else
        {
            RatesOfChange(k1, x, t);
            for (int i = 0; i < x.Length; i++)
            {
                P[i] = x[i] + (h / 24.0) *
                    (55.0 * k1[i] - 59.0 * dm1[i] + 37.0 * dm2[i] - 9.0 * dm3[i]);
            }
            RatesOfChange(dp1, P, t + h);
            abmRms2 = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                store[i] = x[i];
                x[i] += (h / 24.0) * (9 * dp1[i] + 19.0 * k1[i] - 5.0 * dm1[i] + dm2[i]);
                dm3[i] = dm2[i];
                dm2[i] = dm1[i];
                dm1[i] = k1[i];
                ym3[i] = ym2[i];
                ym2[i] = ym1[i];
                ym1[i] = store[i];
                abmRms2 += (x[i] - P[i]) * (x[i] - P[i]) / (x[i] + P[i]) / (x[i] + P[i]);
            }
            abmRms2 /= x.Length;
            if (abmSteps < 5) abmSteps += 1;
            return t + h;
        }
    }

    public double abmError()
    {
        return abmRms2;
    }
}
