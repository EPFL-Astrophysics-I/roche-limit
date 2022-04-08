using UnityEngine;

public class CircularOrbitSlideController : SimulationSlideController
{
    private CircularOrbitSimulation sim;
    private CircularOrbitPrefabs prefabs;

    [SerializeField] private bool useLights = true;
    [SerializeField] private bool rocheLimit = true;

    private void Awake()
    {
        sim = (CircularOrbitSimulation)simulation;
        if (!sim.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Could not find CircularOrbitPrefabs component");
        }
    }

    public override void InitializeSlide()
    {
        if (!prefabs)
        {
            return;
        }

        prefabs.SetLightsVisibility(useLights);
        prefabs.SetRocheLimitVisibility(rocheLimit);
    }
}
