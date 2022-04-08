using UnityEngine;

[RequireComponent(typeof(TidalDistortionPrefabs))]
public class TidalDistortionSimulation : Simulation
{
    private TidalDistortionPrefabs prefabs;

    [Header("Units")]
    [SerializeField] private Units.UnitLength unitLength = Units.UnitLength.EarthRadius;
    [SerializeField] private Units.UnitMass unitMass = Units.UnitMass.EarthMass;
    [SerializeField] private Units.UnitTime unitTime = Units.UnitTime.Day;
    [SerializeField, Min(0)] private float timeScale = 1;
    private double newtonG;

    [Header("Parameters")]
    [SerializeField] private float primaryRadius = 3f;
    [SerializeField] private float primaryDensity = 1f;
    [SerializeField] private float satelliteDistance = 10f;
    [SerializeField] private float satelliteRadius = 0.5f;
    [SerializeField] private float satelliteDensity = 1f;
    private double primaryMass;

    [Header("Solver")]
    [SerializeField, Min(1)] private int numSubsteps = 1;

    // Constant
    private float r;
    private float omega;

    // Variable
    private float theta;

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
            prefabs.InstantiateAllPrefabs(primaryRadius, primaryDensity, satelliteDistance, satelliteRadius, satelliteDensity);
        }
    }

    private void Start()
    {
        // Compute initial conditions
        float rx = prefabs.satellite.position.x;
        float rz = prefabs.satellite.position.z;
        r = Mathf.Sqrt(rx * rx + rz * rz);
        theta = Mathf.Atan2(rz, rx);
        omega = Mathf.Sqrt((float)(newtonG * primaryMass / r / r / r));
    }

    private void FixedUpdate()
    {
        if (timeScale == 0)
        {
            return;
        }

        float substep = timeScale * Time.fixedDeltaTime / numSubsteps;
        for (int i = 0; i < numSubsteps; i++)
        {
            StepForward(substep);
        }

        // Update game object positions
        float rx = r * Mathf.Cos(theta);
        float rz = r * Mathf.Sin(theta);
        prefabs.satellite.position = new Vector3(rx, 0, rz);
    }

    private void StepForward(float deltaTime)
    {
        theta += omega * deltaTime;
        if (theta >= 2 * Mathf.PI)
        {
            theta -= 2 * Mathf.PI;
        }
    }
}
