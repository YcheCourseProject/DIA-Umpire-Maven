/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.PSMDataStructure;

import java.io.Serializable;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PeptideFragment implements Serializable{
    private static final long serialVersionUID = 204820482047L;

    public String IonType="";
    public float ObservedMZ;
    public float FragMZ;
    public int Charge;

    public String GetFragKey() {
        return IonType + "_" + Charge;
    }
}
