using UnityEngine;

[RequireComponent(typeof(CircularOrbitPrefabs))]
public class CircularOrbitSimulation : Simulation
{
    private CircularOrbitPrefabs prefabs;

    [Header("Units")]
    [SerializeField] private Units.UnitLength unitLength = Units.UnitLength.EarthRadius;
    [SerializeField] private Units.UnitMass unitMass = Units.UnitMass.EarthMass;
    [SerializeField] private Units.UnitTime unitTime = Units.UnitTime.Day;
    [SerializeField, Min(0)] private float timeScale = 1;
    private double newtonG;

    [Header("Parameters")]
    [SerializeField] private int numParticles = 100;
    [SerializeField] private float primaryRadius = 3f;
    [SerializeField] private float ringRadius = 4;
    [SerializeField] private float ringDispersion = 0.5f;
    [SerializeField] private float particleRadius = 0.1f;
    [SerializeField] private float satelliteDistance = 10;
    [SerializeField] private float satelliteRadius = 0.5f;
    [SerializeField] private float rocheRadius = 6;
    private double primaryMass;
    private int numBodies;

    [Header("Solver")]
    [SerializeField, Min(1)] private int numSubsteps = 1;

    // Constant
    private float[] r;
    private float[] omega;

    // Variable
    private float[] theta;

    private void Awake()
    {
        // Set the gravitational constant for these units
        newtonG = Units.NewtonG(unitTime, unitLength, unitMass);
        //Debug.Log("G = " + G);
        primaryMass = Units.EarthMass(unitMass);
        //Debug.Log("M = " + M);

        // Create all objects with assigned prefabs
        if (transform.TryGetComponent(out prefabs))
        {
            prefabs.InstantiateAllPrefabs(numParticles, primaryRadius, ringRadius, ringDispersion,
                particleRadius, satelliteDistance, satelliteRadius, rocheRadius);
        }
    }

    private void Start()
    {
        numBodies = prefabs.ring.childCount;
        if (prefabs.satellite)
        {
            numBodies += 1;
        }

        r = new float[numBodies];
        omega = new float[numBodies];
        theta = new float[numBodies];

        // Compute initial conditions
        int startIndex = 0;
        if (prefabs.satellite)
        {
            float rx = prefabs.satellite.position.x;
            float rz = prefabs.satellite.position.z;
            r[0] = Mathf.Sqrt(rx * rx + rz * rz);
            theta[0] = Mathf.Atan2(rz, rx);
            omega[0] = Mathf.Sqrt((float)(newtonG * primaryMass / r[0] / r[0] / r[0]));
            startIndex = 1;
        }

        for (int i = startIndex; i < numBodies; i++)
        {
            Transform particle = prefabs.ring.GetChild(i - startIndex);
            float rx = particle.position.x;
            float rz = particle.position.z;
            r[i] = Mathf.Sqrt(rx * rx + rz * rz);
            theta[i] = Mathf.Atan2(rz, rx);
            omega[i] = Mathf.Sqrt((float)(newtonG * primaryMass / r[i] / r[i] / r[i]));
        }
    }

    private void FixedUpdate()
    {
        float substep = timeScale * Time.fixedDeltaTime / numSubsteps;
        for (int i = 0; i < numSubsteps; i++)
        {
            StepForward(substep);
        }

        // Update game object positions
        int startIndex = 0;
        if (prefabs.satellite)
        {
            float rx = r[0] * Mathf.Cos(theta[0]);
            float rz = r[0] * Mathf.Sin(theta[0]);
            prefabs.satellite.position = new Vector3(rx, 0, rz);
            startIndex = 1;
        }

        for (int i = startIndex; i < numBodies; i++)
        {
            Transform particle = prefabs.ring.GetChild(i - startIndex);
            float rx = r[i] * Mathf.Cos(theta[i]);
            float rz = r[i] * Mathf.Sin(theta[i]);
            particle.position = new Vector3(rx, 0, rz);
        }
    }

    private void StepForward(float deltaTime)
    {
        for (int i = 0; i < numBodies; i++)
        {
            theta[i] += omega[i] * deltaTime;
            if (theta[i] >= 2 * Mathf.PI)
            {
                theta[i] -= 2 * Mathf.PI;
            }
        }
    }
}
