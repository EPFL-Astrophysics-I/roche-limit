using UnityEngine;

public class RocheBandSlideController : SimulationSlideController
{
    private RocheBandSimulation sim;
    private RocheBandPrefabs prefabs;

    [SerializeField] private bool useLights = true;

    private void Awake()
    {
        sim = (RocheBandSimulation)simulation;
        if (!sim.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Could not find RocheBandPrefabs component");
        }
    }

    public override void InitializeSlide()
    {
        if (!prefabs)
        {
            return;
        }

        prefabs.SetLightsVisibility(useLights);
        SetRigidLimitVisibility(true);
        SetFluidLimitVisibility(false);
    }

    public void SetRigidLimitVisibility(bool visible)
    {
        if (prefabs.rigidLimitLR)
        {
            prefabs.rigidLimitLR.gameObject.SetActive(visible);
        }
    }

    public void SetFluidLimitVisibility(bool visible)
    {
        if (prefabs.fluidLimitLR)
        {
            prefabs.fluidLimitLR.gameObject.SetActive(visible);
        }
    }
}
