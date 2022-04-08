using UnityEngine;

[RequireComponent(typeof(RocheBandPrefabs))]
public class RocheBandSimulation : Simulation
{
    private RocheBandPrefabs prefabs;

    [Header("Parameters")]
    [SerializeField] private float primaryRadius = 3f;
    [SerializeField] private float densityRatio = 1.5f;
    //[SerializeField] private float rigidLimitRadius = 5;
    //[SerializeField] private float fluidLimitRadius = 7;

    private void Awake()
    {
        // Create all objects with assigned prefabs
        if (transform.TryGetComponent(out prefabs))
        {
            prefabs.InstantiateAllPrefabs(primaryRadius, densityRatio);
        }
    }
}
