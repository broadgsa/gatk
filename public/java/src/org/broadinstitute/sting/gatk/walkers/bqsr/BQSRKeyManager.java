package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.*;

/**
 * This class provides all the functionality for the BitSet representation of the keys to the hash table of BQSR
 *
 * It also handles the event type "covariate" which is not exactly a covariate, but is added as a key to the hashmap. The Key Manager will
 * add the event type as a bitset to the end of the covariate bitset key. This way, it won't get int the way of masking the information
 * out of the key for the actual covariates, and having the covariates handle it. The key manager handles the event type.
 *
 * The keys represented by this key manager will always have the same order:
 *
 * RequiredCovariate1, RequiredCovariate2, ..., RequiredCovariateN, OptionalCovariate1, OptionalCovariateID, EventType
 * RequiredCovariate1, RequiredCovariate2, ..., RequiredCovariateN, OptionalCovariate2, OptionalCovariateID, EventType
 * ...
 * RequiredCovariate1, RequiredCovariate2, ..., RequiredCovariateN, OptionalCovariateN, OptionalCovariateID, EventType
 *
 *
 * Note that Optional Covariates are optional, and the Key Manager should operate without them if necessary.
 *
 * @author Mauricio Carneiro
 * @since 3/6/12
 */
public class BQSRKeyManager {

    private final Covariate[] requiredCovariates;
    private final Covariate[] optionalCovariates;
    private final RequiredCovariateInfo[] requiredCovariatesInfo;
    private final OptionalCovariateInfo[] optionalCovariatesInfo;
    private final Map<String, Short> covariateNameToIDMap;

    private int nRequiredBits;                                                                                          // Number of bits used to represent the required covariates

    private final int optionalCovariateOffset;
    private final int optionalCovariateIDOffset;

    private final long optionalCovariateMask;                                                                           // Standard mask for optional covariates key
    private final long optionalCovariateIDMask;                                                                         // Standard mask for optional covariates order key
    private final long eventIDMask;                                                                                     // Standard mask for event ID

    /**
     * Initializes the KeyManager with the total number of covariates to use
     *
     * @param requiredCovariates the ordered list of required covariates
     * @param optionalCovariates the ordered list of optional covariates
     */
    public BQSRKeyManager(final List<Covariate> requiredCovariates, final List<Covariate> optionalCovariates) {
        this.requiredCovariates = new Covariate[requiredCovariates.size()];
        this.optionalCovariates = new Covariate[optionalCovariates.size()];
        requiredCovariatesInfo = new RequiredCovariateInfo[requiredCovariates.size()];                                  // initialize the required covariates list
        optionalCovariatesInfo = new OptionalCovariateInfo[optionalCovariates.size()];                                  // initialize the optional covariates list (size may be 0, it's okay)
        covariateNameToIDMap = new HashMap<String, Short>(optionalCovariates.size()*2);                                 // the map from covariate name to covariate id (when reading GATK Reports, we get the IDs as names of covariates)
        
        nRequiredBits = 0;
        for (int i = 0; i < requiredCovariates.size(); i++) {                                                           // create a list of required covariates with the extra information for key management
            final Covariate required = requiredCovariates.get(i);
            final int nBits = required.numberOfBits();                                                                  // number of bits used by this covariate
            final long mask = genericMask(nRequiredBits, nBits);                                                        // create a mask for this covariate
            this.requiredCovariates[i] = required;
            requiredCovariatesInfo[i] = new RequiredCovariateInfo(nBits, nRequiredBits, mask, required);                // Create an object for this required covariate
            nRequiredBits += nBits;
        }

        final int bitsInEventType = numberOfBitsToRepresent(EventType.values().length);
        eventIDMask = genericMask(nRequiredBits, bitsInEventType);

        short id = 0;
        int nOptionalBits = 0;
        for (int i = 0; i < optionalCovariates.size(); i++) {
            final Covariate optional = optionalCovariates.get(i);
            nOptionalBits = Math.max(nOptionalBits, optional.numberOfBits());                                           // optional covariates are represented by the number of bits needed by biggest covariate
            this.optionalCovariates[i] = optional;
            optionalCovariatesInfo[i] = new OptionalCovariateInfo(id, optional);
            final String covariateName = optional.getClass().getSimpleName().split("Covariate")[0];                     // get the name of the covariate (without the "covariate" part of it) so we can match with the GATKReport
            covariateNameToIDMap.put(covariateName, id);
            id++;
        }

        optionalCovariateOffset = nRequiredBits + bitsInEventType;
        optionalCovariateMask = genericMask(optionalCovariateOffset, nOptionalBits);                                    // the generic mask to extract optional covariate bits from the combined bitset
        optionalCovariateIDOffset = nRequiredBits + bitsInEventType + nOptionalBits;
        final int nOptionalIDBits = numberOfBitsToRepresent(optionalCovariates.size());                                 // number of bits used to represent the covariate ID
        optionalCovariateIDMask = genericMask(optionalCovariateIDOffset, nOptionalIDBits);                              // the generic mask to extract optional covariate ID bits from the combined bitset

        final int totalNumberOfBits = optionalCovariateIDOffset + nOptionalIDBits;                                      // total number of bits used in the final key
        if ( totalNumberOfBits > 64 )
            throw new UserException.BadInput("The total number of bits used for the master BQSR key is greater than 64 and cannot be represented in a Long");
    }

    /**
     * Generates one key per optional covariate.
     * 
     * Keys include all required covariates, the standard covariate and the event type.
     *
     * Example allKeys:
     * RG, QUAL, CYCLE, CONTEXT
     *
     * List of BitSets returned by this example (given eventType):
     * RG, QUAL, CYCLE, EVENT
     * RG, QUAL, CONTEXT, EVENT
     *
     * Note: If there are no optional covariates, only one bitset key will be returned with all the required covariates and the event type
     *
     * @param allKeys      The keys in long representation for each covariate
     * @param eventType The type of event described by this keyset (e.g. mismatches, insertions, deletions)
     * @return one key in long representation per covariate
     */
    public Long[] longsFromAllKeys(final Long[] allKeys, final EventType eventType) {
        final ArrayList<Long> allFinalKeys = new ArrayList<Long>();                                                     // Generate one key per optional covariate

        int covariateIndex = 0;
        long masterKey = 0L;                                                                                            // This will be a master key holding all the required keys, to replicate later on
        for (RequiredCovariateInfo infoRequired : requiredCovariatesInfo)
            masterKey |= (allKeys[covariateIndex++] << infoRequired.offset);

        final long eventKey = keyFromEvent(eventType);                                                                  // create a key for the event type
        masterKey |= (eventKey << nRequiredBits);

        for (int i = 0; i < optionalCovariatesInfo.length; i++) {
            final Long covariateKey = allKeys[covariateIndex++];
            if (covariateKey == null)
                continue;                                                                                               // do not add nulls to the final set of keys.

            long newKey = masterKey | (covariateKey << optionalCovariateOffset);
            newKey |= (optionalCovariatesInfo[i].covariateID << optionalCovariateIDOffset);

            allFinalKeys.add(newKey);                                                                                   // add this key to the list of keys
        }

        if (optionalCovariatesInfo.length == 0)                                                                         // special case when we have no optional covariates
            allFinalKeys.add(masterKey);

        return allFinalKeys.toArray(new Long[allFinalKeys.size()]);
    }

    /**
     * Generates one key for the covariates represented in Object[] key
     *
     * The covariates will have the actual objects produced by the covariates (probably read from the recalibration data file)
     * and will contain all required covariates and one (or none) optional covariates. Therefore, the product is one key, not many.
     *
     * Example key:
     * RG, QUAL, CYCLE, CYCLE_ID, EventType
     *
     * @param key list of objects produced by the required covariates followed by one or zero optional covariates.
     * @return a key representing these objects.
     */
    public Long longFromKey(Object[] key) {
        int requiredCovariate = 0;
        long masterKey = 0L;                                                                                            // This will be a master key holding all the required keys, to replicate later on
        for (RequiredCovariateInfo infoRequired : requiredCovariatesInfo)
            masterKey |= (infoRequired.covariate.longFromKey(key[requiredCovariate++]) << infoRequired.offset);

        final int eventIndex = key.length - 1;                                                                          // the event type is always the last key
        final long eventKey = keyFromEvent((EventType) key[eventIndex]);                                                // create a key for the event type
        masterKey |= (eventKey << nRequiredBits);

        if (optionalCovariatesInfo.length > 0) {
            final int covariateIndex = requiredCovariatesInfo.length;                                                   // the optional covariate index in the key array
            final int covariateIDIndex = covariateIndex + 1;                                                            // the optional covariate ID index is right after the optional covariate's
            final short covariateID = parseCovariateID(key[covariateIDIndex]);                                          // when reading the GATK Report the ID may come in a String instead of an index
            final OptionalCovariateInfo infoOptional = optionalCovariatesInfo[covariateID];                             // so we can get the optional covariate information

            final long covariateKey = infoOptional.covariate.longFromKey(key[covariateIndex]);                          // convert the optional covariate key into a bitset using the covariate's interface
            masterKey |= (covariateKey << optionalCovariateOffset);
            masterKey |= (infoOptional.covariateID << optionalCovariateIDOffset);
        }

        return masterKey;
    }

    /**
     * Covariate id can be either the covariate name (String) or the actual id (short). This method
     * finds it's type and converts accordingly to the short notation.
     *
     * @param id the string or short representation of the optional covariate id
     * @return the short representation of the optional covariate id.
     */
    private short parseCovariateID(final Object id) {
        return (id instanceof String) ? covariateNameToIDMap.get(id.toString()) : (Short) id;
    }

    /**
     * Generates a key set of objects from a combined master key.
     *
     * Masks out each covariate independently and decodes their values (Object) into a keyset
     *
     * @param master the master representation of the keys
     * @return an object array with the values for each key
     */
    public List<Object> keySetFrom(final Long master) {
        final List<Object> objectKeys = new ArrayList<Object>();
        for (RequiredCovariateInfo info : requiredCovariatesInfo) {
            final Long covariateKey = extractKeyFromMaster(master, info.mask, info.offset);                             // get the covariate's key
            objectKeys.add(info.covariate.formatKey(covariateKey));                                                     // convert the key to object using covariate's interface
        }

        if (optionalCovariatesInfo.length > 0) {
            final Long covKey = extractKeyFromMaster(master, optionalCovariateMask, optionalCovariateOffset);           // get the covariate's key
            final int covIDKey = (int)extractKeyFromMaster(master, optionalCovariateIDMask, optionalCovariateIDOffset); // get the covariate's id (to identify which covariate this is)
            Covariate covariate = optionalCovariatesInfo[(short)covIDKey].covariate;                                    // get the corresponding optional covariate object
            objectKeys.add(covariate.formatKey(covKey));                                                                // add the optional covariate key to the key set
            objectKeys.add(covariate.getClass().getSimpleName().split("Covariate")[0]);                                 // add the covariate name using the id
        }

        objectKeys.add(EventType.eventFrom((int)extractKeyFromMaster(master, eventIDMask, nRequiredBits)));             // add the event type object to the key set

        return objectKeys;
    }

    public Covariate[] getRequiredCovariates() {
        return requiredCovariates;
    }

    public Covariate[] getOptionalCovariates() {
        return optionalCovariates;
    }

    public int getNumRequiredCovariates() {
        return requiredCovariates.length;
    }

    public int getNumOptionalCovariates() {
        return optionalCovariates.length;
    }

    /**
     * Creates a mask for the requested covariate to extract the relevant key from a combined master key
     *
     * @param offset  the offset into the master key
     * @param nBits   the number of bits needed by the Covariate to represent its values
     * @return the mask relevant to the covariate
     */
    private long genericMask(final int offset, final int nBits) {
        long mask = 0L;
        for ( int i = 0; i < nBits; i++ )
            mask |= 1L << (offset+i);
        return mask;
    }

    private long extractKeyFromMaster(final long master, final long mask, final int offset) {
        long key = master & mask;
        return key >> offset;
    }

    // cache the key representing an event since it's otherwise created a massive amount of times
    private static final Map<EventType, Long> eventTypeCache = new HashMap<EventType, Long>(EventType.values().length); // event IDs must be longs so that bit-fiddling works
    static {
        for (final EventType eventType : EventType.values())
            eventTypeCache.put(eventType, (long)eventType.index);
    }

    private Long keyFromEvent(final EventType eventType) {
        return eventTypeCache.get(eventType);
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof BQSRKeyManager))
            return false;

        BQSRKeyManager other = (BQSRKeyManager) o;
        if (this == other)
            return true;

        if (requiredCovariatesInfo.length != other.requiredCovariatesInfo.length ||
                optionalCovariatesInfo.length != other.optionalCovariatesInfo.length)
            return false;

        for (int i = 0; i < requiredCovariates.length; i++) {
            Covariate myRequiredCovariate = requiredCovariates[i];
            Covariate otherRequiredCovariate = other.requiredCovariates[i];
            String thisName = myRequiredCovariate.getClass().getSimpleName();
            String otherName = otherRequiredCovariate.getClass().getSimpleName();
            if (!thisName.equals(otherName))
                return false;
        }

        for (int i = 0; i < optionalCovariates.length; i++) {
            Covariate myOptionalCovariate = optionalCovariates[i];
            Covariate otherOptionalCovariate = other.optionalCovariates[i];
            String thisName = myOptionalCovariate.getClass().getSimpleName();
            String otherName = otherOptionalCovariate.getClass().getSimpleName();
            if (!thisName.equals(otherName))
                return false;
        }

        return true;
    }

    /**
     * Calculates the number of bits necessary to represent a given number of elements
     *
     * @param numberOfElements the number of elements to represent (must be positive)
     * @return the number of bits necessary to represent this many elements
     */
    public static int numberOfBitsToRepresent(long numberOfElements) {
        if (numberOfElements < 0)
            throw new ReviewedStingException("Number of elements must be positive: " + numberOfElements);

        if (numberOfElements == 1L)
            return 1;   // special case

        int n = 0;
        numberOfElements--;
        while (numberOfElements > 0) {
            numberOfElements = numberOfElements >> 1;
            n++;
        }
        return n;
    }

    /**
     * Aggregate information for each Covariate
     */
    private static class RequiredCovariateInfo {
        public final int nBits;                                                                                         // number of bits for this key
        public final int offset;                                                                                        // the offset into the master key
        public final long mask;                                                                                         // the mask to pull out this covariate from the combined bitset key ( a mask made from bitsBefore and nBits )
        public final Covariate covariate;                                                                               // this allows reverse lookup of the Covariates in order

        RequiredCovariateInfo(final int nBits, final int offset, final long mask, final Covariate covariate) {
            this.nBits = nBits;
            this.offset = offset;
            this.mask = mask;
            this.covariate = covariate;
        }
    }

    private static class OptionalCovariateInfo {
        public final long covariateID;                                                                                  // cache the covariate ID (must be a long so that bit-fiddling works)
        public final Covariate covariate;

        OptionalCovariateInfo(final long covariateID, final Covariate covariate) {
            this.covariateID = covariateID;
            this.covariate = covariate;
        }
    }
    
}
