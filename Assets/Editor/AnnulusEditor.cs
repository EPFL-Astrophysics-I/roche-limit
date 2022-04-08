using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(Annulus))]
public class AnnulusEditor : Editor
{
    public override void OnInspectorGUI()
    {
        Annulus annulus = (Annulus)target;

        if (DrawDefaultInspector())
        {
            if (annulus.autoUpdate && !Application.isPlaying)
            {
                annulus.Generate();
            }
        }

        GUILayout.BeginHorizontal();
        if (GUILayout.Button("Generate"))
        {
            if (!Application.isPlaying)
            {
                annulus.Generate();
            }
        }
        if (GUILayout.Button("Clear"))
        {
            if (!Application.isPlaying)
            {
                annulus.Clear();
            }
        }
        GUILayout.EndHorizontal();
    }
}
