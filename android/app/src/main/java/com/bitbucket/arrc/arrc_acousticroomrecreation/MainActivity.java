package com.bitbucket.arrc.arrc_acousticroomrecreation;

import android.app.Activity;
import android.content.Intent;
import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.Toast;

import org.json.JSONException;
import org.json.JSONObject;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.net.HttpURLConnection;
import java.net.URL;

public class MainActivity extends Activity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        Button goButton = findViewById(R.id.main_button_go);
        Button dokmanicButton = findViewById(R.id.main_button_dokmanic);

        final EditText distSpkrET = findViewById(R.id.main_speaker_dist);
        final EditText heightDiffET = findViewById(R.id.main_height_diff);
        final EditText ipET = findViewById(R.id.main_ip);

        goButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                Bundle bundle = new Bundle();
                try {
                    bundle.putFloat("distSpkr", Float.valueOf(distSpkrET.getText().toString()));
                    bundle.putFloat("heightDiff", Float.valueOf(heightDiffET.getText().toString()));
                    String ip = ipET.getText().toString();
                    if (ip.equals("")) {
                        throw new IOException();
                    } else {
                        bundle.putString("ip", ip);
                    }
                } catch (NumberFormatException|IOException e) {
                    Toast.makeText(v.getContext(),
                            "Please enter valid numbers and IP.", Toast.LENGTH_SHORT).show();
                    return;
                }
                Intent intent = new Intent(v.getContext(), SweepActivity.class);
                intent.putExtras(bundle);
                startActivity(intent);
            }
        });

        dokmanicButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                Bundle bundle = new Bundle();
                int[] pk1 = {1102, 1343, 1479, 1643, 1798, 1927, 2070, 2131, 2249, 2616, 2734, 2958, 3038, 3139, 3212, 3390, 3495, 3974};
                bundle.putIntArray("peaks1", pk1);
                int[] pk2 = {1168, 1403, 1510, 1683, 1844, 2069, 2212, 2240, 2370, 2545, 2808, 2898, 3106, 3155, 3292, 3355, 3438, 3618, 3970};
                bundle.putIntArray("peaks2", pk2);
                int[] pk3 = {1111, 1376, 1451, 1632, 1808, 2006, 2161, 2227, 2317, 2358, 2653, 2709, 2847, 3095, 3196, 3284, 3386, 3998};
                bundle.putIntArray("peaks3", pk3);
                int[] pk4 = {1048, 1299, 1445, 1566, 1730, 1975, 2116, 2279, 2376, 2748, 2834, 2916, 3027, 3144, 3337, 3496, 3868, 3913};
                bundle.putIntArray("peaks4", pk4);
                int[] pk5 = {1083, 1338, 1454, 1616, 1781, 1948, 2097, 2174, 2299, 2377, 2648, 2740, 2857, 2907, 3109, 3178, 3350, 3445, 3889, 3948};
                bundle.putIntArray("peaks5", pk5);
                bundle.putInt("fs", 96000);
                String ip = ipET.getText().toString();
                if (ip.equals("")) {
                    Toast.makeText(v.getContext(),
                            "Please enter valid numbers and IP.", Toast.LENGTH_SHORT).show();
                } else {
                    bundle.putString("ip", ip);
                    Intent intent = new Intent(v.getContext(), RoomShapeActivity.class);
                    intent.putExtras(bundle);
                    startActivity(intent);
                }
            }
        });
    }
}
