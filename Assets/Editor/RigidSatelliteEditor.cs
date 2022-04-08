using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(RigidSatellite))]
public class RigidSatelliteEditor : Editor
{
    public override void OnInspectorGUI()
    {
        RigidSatellite rigidSatellite = (RigidSatellite)target;

        if (DrawDefaultInspector())
        {
            if (rigidSatellite.autoUpdate && !Application.isPlaying)
            {
                rigidSatellite.Generate();
            }
        }

        GUILayout.BeginHorizontal();
        if (GUILayout.Button("Generate"))
        {
            if (!Application.isPlaying)
            {
                rigidSatellite.Generate();
            }
        }
        if (GUILayout.Button("Clear"))
        {
            if (!Application.isPlaying)
            {
                rigidSatellite.Clear();
            }
        }
        GUILayout.EndHorizontal();
    }
}
