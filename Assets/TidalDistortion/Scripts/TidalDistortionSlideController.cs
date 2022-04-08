using UnityEngine;

public class TidalDistortionSlideController : SimulationSlideController
{
    private TidalDistortionSimulation sim;
    private TidalDistortionPrefabs prefabs;

    [SerializeField] private bool useLights = true;
    [SerializeField] private bool rocheLimit = true;

    private void Awake()
    {
        sim = (TidalDistortionSimulation)simulation;
        if (!sim.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Could not find TidalDistortionPrefabs component");
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
